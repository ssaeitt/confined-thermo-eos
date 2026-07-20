classdef FlashEngine < handle
    % FLASHENGINE Multi-variable non-isobaric VLE flash calculation engine.
    % Traces dew-point phase boundaries under nanoporous confinement using
    % hybrid Newton-Raphson and Broyden Quasi-Newton solvers with Armijo
    % backtracking line search and Euclidean domain projections.

    properties (SetAccess = private)
        EOS thermo.ConfinedEOS         % Handle to active thermodynamic operator
        Stability solvers.StabilityTester % Handle to stability analysis engine

        % Numerical Solver Configuration Controls
        MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 1000
        Tolerance (1,1) double {mustBePositive, mustBeReal} = 1e-8
        JacEpsilon (1,1) double {mustBePositive, mustBeReal} = 1e-6
        LineSearchC1 (1,1) double {mustBePositive, mustBeReal} = 1e-4
        MaxLineSearch (1,1) double {mustBeInteger, mustBePositive} = 10
        JacRegularizer (1,1) double {mustBePositive, mustBeReal} = 1e-8
    end

    methods
        function obj = FlashEngine(eosEngine, stabilityEngine, varargin)
            % Constructor attaching thermodynamic and stability operators
            validateattributes(eosEngine, {'thermo.ConfinedEOS'}, {'scalar'});
            validateattributes(stabilityEngine, {'solvers.StabilityTester'}, {'scalar'});

            obj.EOS = eosEngine;
            obj.Stability = stabilityEngine;

            if nargin > 2
                obj.parseSolverOptions(varargin{:});
            end
        end

        function [Pdew, K_final, Pcap_final, Pliq_final, solverStats] = solveDewPoint(obj, T, P_guess, z, r_cap, varargin)
            % Main execution entry point for non-isobaric dew point prediction.
            % Automatically handles mode switching between confined ('Pc') and bulk ('PL') states.

            z = reshape(z, [], 1);
            z = z / sum(z);

            % Parse optional configuration flags for this run
            solverType = 'newton';     % Default to full Newton-Raphson
            capMode    = 'Pc';         % Default to Capillary Pressure tracking
            K_seed     = [];
            Pcap_seed  = [];

            for idx = 1:2:length(varargin)
                switch lower(string(varargin{idx}))
                    case "solver",    solverType = lower(varargin{idx+1});
                    case "capmode",   capMode    = varargin{idx+1};
                    case "k_seed",    K_seed     = varargin{idx+1};
                    case "pcap_seed", Pcap_seed  = varargin{idx+1};
                end
            end

            % 1. Automatic Regime Detection & Singularity Guard
            if isinf(r_cap)
                capMode = 'PL'; % Force liquid pressure tracking in unconfined bulk media
            end

            % 2. Initialize Thermodynamic State Vector Seeding
            if isempty(K_seed) || isempty(Pcap_seed)
                [isUnstable, w_trial, ~, K_tpd, ~, PL_tpd, Pc_tpd] = ...
                    obj.Stability.evaluateDewStability(P_guess, T, z, r_cap);

                if isUnstable && ~any(isnan(K_tpd))
                    K_init = K_tpd;
                    Pcap_init = Pc_tpd;
                    PL_init   = PL_tpd;
                else
                    % Fallback to Wilson empirical seeding if feed appears locally stable
                    Pc_vec = obj.EOS.Fluid.Pc(:);
                    Tc_vec = obj.EOS.Fluid.Tc(:);
                    om_vec = obj.EOS.Fluid.omega(:);
                    K_init = (Pc_vec ./ P_guess) .* exp(5.37 .* (1 + om_vec) .* (1 - Tc_vec ./ T));
                    Pcap_init = 0.0;
                    PL_init   = P_guess;
                end
            else
                K_init    = K_seed;
                Pcap_init = Pcap_seed;
                PL_init   = max(P_guess - Pcap_seed, 1e3);
            end

            % Assemble log-explicit independent variable vector u
            lnK0  = log(max(K_init(:), 1e-14));
            lnPV0 = log(max(P_guess, 1e5));

            if strcmpi(capMode, 'pl')
                lnPmech0 = log(max(PL_init, 1e5));
            else
                lnPmech0 = log(max(Pcap_init, 10.0)); % Clamp minimum Pc seed to 10 Pa
            end

            u0 = [lnK0; lnPV0; lnPmech0];

            % Establish local search boundaries matching legacy robust frameworks
            lnPV_min = log(0.5 * P_guess);
            lnPV_max = log(3.0 * P_guess);

            % 3. Execute Selected Iterative Solver Scheme
            if strcmpi(solverType, 'quasinewton')
                [u_converged, stats] = obj.executeBroydenSolver(u0, T, z, r_cap, capMode, lnPV_min, lnPV_max);
            else
                [u_converged, stats] = obj.executeNewtonSolver(u0, T, z, r_cap, capMode, lnPV_min, lnPV_max);
            end

            % 4. Unpack Converged State Boundaries
            nc = obj.EOS.Fluid.NC;
            K_final  = exp(u_converged(1:nc));
            Pdew     = exp(u_converged(nc+1));

            if strcmpi(capMode, 'pl')
                Pliq_final = exp(u_converged(nc+2));
                Pcap_final = max(Pdew - Pliq_final, 0.0);
            else
                Pcap_final = exp(u_converged(nc+2));
                Pliq_final = max(Pdew - Pcap_final, 1e3);
            end

            solverStats = stats;
        end
    end

    methods (Access = private)
        function [u, stats] = executeNewtonSolver(obj, u0, T, z, r_cap, mode, lnPV_min, lnPV_max)
            % Full Newton-Raphson algorithm with numerical Jacobian and backtracking line search
            u = obj.projectFeasibleDomain(u0, mode, r_cap, lnPV_min, lnPV_max);
            iterLog = zeros(obj.MaxIterations, 1);

            for iter = 1:obj.MaxIterations
                % Evaluate baseline residual
                [res, ~] = obj.evaluateResidualVector(u, T, z, r_cap, mode);
                norm_res = norm(res, inf);
                iterLog(iter) = norm_res;

                if norm_res < obj.Tolerance
                    stats = struct('iterations', iter, 'residual_norm', norm_res, 'converged', true, 'history', iterLog(1:iter));
                    return;
                end

                % Compute Finite-Difference Jacobian Matrix
                J = obj.computeProjectedCentralJacobian(u, T, z, r_cap, mode, lnPV_min, lnPV_max);

                if obj.JacRegularizer > 0
                    J = J + obj.JacRegularizer * eye(size(J));
                end

                step = -J \ res;

                % Armijo Backtracking Line Search
                alpha = 1.0;
                u_new = obj.projectFeasibleDomain(u + alpha * step, mode, r_cap, lnPV_min, lnPV_max);
                [res_new, ~] = obj.evaluateResidualVector(u_new, T, z, r_cap, mode);

                ls_iter = 0;
                while norm(res_new, inf) > (1 - obj.LineSearchC1 * alpha) * norm_res && ls_iter < obj.MaxLineSearch
                    alpha = alpha * 0.5;
                    u_new = obj.projectFeasibleDomain(u + alpha * step, mode, r_cap, lnPV_min, lnPV_max);
                    [res_new, ~] = obj.evaluateResidualVector(u_new, T, z, r_cap, mode);
                    ls_iter = ls_iter + 1;
                end

                % Check step size stalling
                if norm(u_new - u, inf) < 1e-12
                    break;
                end
                u = u_new;
            end

            stats = struct('iterations', iter, 'residual_norm', norm_res, 'converged', false, 'history', iterLog(1:iter));
        end

        function [u, stats] = executeBroydenSolver(obj, u0, T, z, r_cap, mode, lnPV_min, lnPV_max)
            % Broyden Quasi-Newton rank-1 update solver
            u = obj.projectFeasibleDomain(u0, mode, r_cap, lnPV_min, lnPV_max);
            [res, ~] = obj.evaluateResidualVector(u, T, z, r_cap, mode);
            norm_res = norm(res, inf);

            % Seed initial Jacobian via finite differences
            J = obj.computeProjectedCentralJacobian(u, T, z, r_cap, mode, lnPV_min, lnPV_max);
            B = J;

            if obj.JacRegularizer > 0
                B = B + obj.JacRegularizer * eye(size(B));
            end
            iterLog = zeros(obj.MaxIterations, 1);

            for iter = 1:obj.MaxIterations
                iterLog(iter) = norm_res;
                if norm_res < obj.Tolerance
                    stats = struct('iterations', iter, 'residual_norm', norm_res, 'converged', true, 'history', iterLog(1:iter));
                    return;
                end

                step = -B \ res;

                % Line search
                alpha = 1.0;
                u_new = obj.projectFeasibleDomain(u + alpha * step, mode, r_cap, lnPV_min, lnPV_max);
                [res_new, ~] = obj.evaluateResidualVector(u_new, T, z, r_cap, mode);

                ls_iter = 0;
                while norm(res_new, inf) > (1 - obj.LineSearchC1 * alpha) * norm_res && ls_iter < obj.MaxLineSearch
                    alpha = alpha * 0.5;
                    u_new = obj.projectFeasibleDomain(u + alpha * step, mode, r_cap, lnPV_min, lnPV_max);
                    [res_new, ~] = obj.evaluateResidualVector(u_new, T, z, r_cap, mode);
                    ls_iter = ls_iter + 1;
                end

                % Broyden Rank-1 Matrix Update: B_{k+1} = B_k + ((y - B_k * s) * s^T) / (s^T * s)
                s = u_new - u;
                y_vec = res_new - res;
                if norm(s) > 1e-14
                    B = B + ((y_vec - B * s) * s') / (s' * s);
                end

                if norm(s, inf) < 1e-12
                    break;
                end

                u = u_new;
                res = res_new;
                norm_res = norm(res, inf);
            end

            stats = struct('iterations', iter, 'residual_norm', norm_res, 'converged', false, 'history', iterLog(1:iter));
        end

        function [F, stateData] = evaluateResidualVector(obj, u, T, z, r_cap, mode)
            % Evaluates the NC+2 residual vector for non-isobaric phase equilibrium
            nc = obj.EOS.Fluid.NC;

            try
                % 1. Unpack state coordinates
                K    = exp(u(1:nc));
                P_v  = exp(u(nc+1));

                if strcmpi(mode, 'pl')
                    P_l   = exp(u(nc+2));
                    P_cap = max(P_v - P_l, 0.0);
                else
                    P_cap = exp(u(nc+2));
                    P_l   = max(P_v - P_cap, 1e3);
                end

                % 2. Calculate trial liquid composition for dew-point tracing (x_i = z_i / K_i)
                x_raw = z ./ max(K, 1e-14);
                x_raw = max(x_raw, 1e-16);
                sx = sum(x_raw);

                if ~isfinite(sx) || sx <= 0
                    F = 1e3 * ones(nc + 2, 1);
                    stateData = struct('P_v', P_v, 'P_l', P_v, 'P_cap', 0, 'Z_V', 1, 'Z_L', 1, 'x_norm', z);
                    return;
                end
                x_norm = x_raw / sx;

                % 3. Evaluate Phase Fugacities and Molar Volumes via ConfinedEOS
                [lnphi_V, Z_V, V_V] = obj.EOS.calculateState(P_v, T, z, -1, r_cap);
                [lnphi_L, Z_L, V_L] = obj.EOS.calculateState(P_l, T, x_norm, 1, r_cap);

                % 4. Evaluate MacLeod-Sugden Parachor Capillary Discontinuity
                rho_l_SI = 1.0 / V_L; % [mol/m^3]
                rho_v_SI = 1.0 / V_V; % [mol/m^3]

                Pcap_pred = obj.evaluateYoungLaplaceIFT(x_norm, z, rho_l_SI, rho_v_SI, r_cap);

                % 5. Assemble Residual Equations
                F = zeros(nc + 2, 1);

                % Chemical potential equality (Iso-fugacity with pressure jump):
                F(1:nc) = log(K) - (lnphi_L - lnphi_V + log(P_l) - log(P_v));

                % Stoichiometric mass balance (Dew point target): sum(z_i / K_i) - 1 = 0
                F(nc+1) = sx - 1.0;

                % Mechanical boundary closure
                if strcmpi(mode, 'pl')
                    F(nc+2) = P_v - P_l - Pcap_pred;
                else
                    F(nc+2) = Pcap_pred - P_cap;
                end

                % Final check to prevent any accidental unphysical NaN leakage
                if any(~isfinite(F))
                    F = 1e3 * ones(nc + 2, 1);
                end
                stateData = struct('P_v', P_v, 'P_l', P_l, 'P_cap', P_cap, 'Z_V', Z_V, 'Z_L', Z_L, 'x_norm', x_norm);

            catch
                % Total fallback isolation envelope to ensure line-search backtrack triggered cleanly
                F = 1e3 * ones(nc + 2, 1);
                stateData = struct('P_v', exp(u(nc+1)), 'P_l', exp(u(nc+1)), 'P_cap', 0, 'Z_V', 1, 'Z_L', 1, 'x_norm', z);
                return;
            end
        end

        function J = computeProjectedCentralJacobian(obj, u, T, z, r_cap, mode, lnPV_min, lnPV_max)
            % Central/Forward finite-difference Jacobian construction
            n = length(u);
            J = zeros(n, n);

            for j = 1:n
                uj = u(j);
                h = obj.JacEpsilon * max(1.0, abs(uj));

                % Formulate perturbed state vectors
                u_f = u; u_f(j) = uj + h;
                u_b = u; u_b(j) = uj - h;

                % Pass perturbations through the feasibility projection block
                u_f_proj = obj.projectFeasibleDomain(u_f, mode, r_cap, lnPV_min, lnPV_max);
                u_b_proj = obj.projectFeasibleDomain(u_b, mode, r_cap, lnPV_min, lnPV_max);

                denom = u_f_proj(j) - u_b_proj(j);

                if abs(denom) < 1e-10
                    % Fall back to a projected forward difference if bounded
                    u_f_proj = obj.projectFeasibleDomain(u, mode, r_cap, lnPV_min, lnPV_max);
                    u_pert = u_f_proj; u_pert(j) = u_f_proj(j) + h;
                    u_pert_proj = obj.projectFeasibleDomain(u_pert, mode, r_cap, lnPV_min, lnPV_max);

                    [res0, ~] = obj.evaluateResidualVector(u_f_proj, T, z, r_cap, mode);
                    [res_pert, ~] = obj.evaluateResidualVector(u_pert_proj, T, z, r_cap, mode);

                    denom_fwd = u_pert_proj(j) - u_f_proj(j);
                    J(:, j) = (res_pert - res0) / max(denom_fwd, 1e-10);
                else
                    [res_f, ~] = obj.evaluateResidualVector(u_f_proj, T, z, r_cap, mode);
                    [res_b, ~] = obj.evaluateResidualVector(u_b_proj, T, z, r_cap, mode);
                    J(:, j) = (res_f - res_b) / denom;
                end
            end
        end

        function Pcap = evaluateYoungLaplaceIFT(obj, x, y, rho_l_SI, rho_v_SI, r_cap)
            % Evaluates interfacial tension and capillary pressure in SI units [Pa]
            if isinf(r_cap)
                Pcap = 0.0;
                return;
            end

            Pch = obj.EOS.Fluid.Parachor(:);

            % Convert molar densities from mol/m^3 to gmol/cm^3 for standard Parachor correlation
            rho_l_cgs = rho_l_SI * 1e-6;
            rho_v_cgs = rho_v_SI * 1e-6;

            ift_param = sum(Pch .* (x .* rho_l_cgs - y .* rho_v_cgs));
            sigma_mN_m = max(ift_param, 0.0)^4;       % Interfacial tension in mN/m
            sigma_N_m  = sigma_mN_m / 1000;            % Convert to N/m

            theta_rad = deg2rad(obj.EOS.Rock.Theta);
            Pcap = (2.0 * sigma_N_m * cos(theta_rad)) / r_cap;
        end

        function u_proj = projectFeasibleDomain(~, u, mode, r_cap, lnPV_min, lnPV_max)
            % Enforces physical boundaries during iteration steps
            u_proj = u;
            n = length(u) - 2;

            % 1. Clamp Vapor Pressure between 0.1 bar (1e4 Pa) and 1000 bar (1e8 Pa)
            u_proj(n+1) = min(max(u_proj(n+1), lnPV_min), lnPV_max);

            % 2. Constrain Mechanical Pressure Variable
            if strcmpi(mode, 'pc')
                margin = log(1 + 1e-5);
                u_proj(n+2) = min(u_proj(n+2), u_proj(n+1) - margin);
            elseif strcmpi(mode, 'pl')
                if isinf(r_cap)
                    % CORRECTED: Binds unconfined liquid pressure boundaries 
                    % to prevent unphysical out-of-bounds parameter execution drift
                    u_proj(n+2) = min(max(u_proj(n+2), lnPV_min), lnPV_max);
                else
                    u_proj(n+2) = min(u_proj(n+2), u_proj(n+1));
                end
            end
        end

        function parseSolverOptions(obj, varargin)
            % Overrides internal numerical tolerances dynamically
            for idx = 1:2:length(varargin)
                switch string(varargin{idx})
                    case "MaxIterations",  obj.MaxIterations  = varargin{idx+1};
                    case "Tolerance",      obj.Tolerance      = varargin{idx+1};
                    case "JacEpsilon",     obj.JacEpsilon     = varargin{idx+1};
                    case "LineSearchC1",   obj.LineSearchC1   = varargin{idx+1};
                    case "JacRegularizer", obj.JacRegularizer = varargin{idx+1};
                end
            end
        end
    end
end