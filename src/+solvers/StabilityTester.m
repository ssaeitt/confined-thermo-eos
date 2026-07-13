classdef StabilityTester < handle
    % STABILITYTESTER Non-isobaric Tangent Plane Distance (TPD) stability engine.
    % Evaluates phase stability under nanoporous confinement by coupling projected
    % gradient descent on the probability simplex with lagged capillary pressure updates.
    
    properties (SetAccess = private)
        EOS thermo.ConfinedEOS  % Handle to active thermodynamic operator
        MaxIterations (1,1) double {mustBeInteger, mustBePositive}
        GradTolerance (1,1) double {mustBePositive, mustBeReal}
        StepTolerance (1,1) double {mustBePositive, mustBeReal}
        Alpha (1,1) double {mustBePositive, mustBeReal}       % Projected gradient step size
        Delta (1,1) double {mustBePositive, mustBeReal}       % Finite-difference perturbation step
        RelaxPL (1,1) double {mustBeReal}                     % Under-relaxation factor for pressure updates
    end
    
    methods
        function obj = StabilityTester(eosEngine, varargin)
            % Constructor attaching physical models and numerical boundaries
            validateattributes(eosEngine, {'thermo.ConfinedEOS'}, {'scalar'});
            obj.EOS = eosEngine;
            
            % Default operational boundaries matching legacy numerics
            obj.MaxIterations = 500;
            obj.GradTolerance = 1e-6;
            obj.StepTolerance = 1e-6;
            obj.Alpha         = 1e-2;
            obj.Delta         = 1e-6;
            obj.RelaxPL       = 0.5;
            
            if nargin > 1
                obj.parseNumericalOverrides(varargin{:});
            end
        end
        
        function [isUnstable, w_star, tpd_min, K_factors, PV, PL, Pc] = evaluateDewStability(obj, Pext, T, z, r_cap)
            % Main entry point evaluating vapor-phase (dew) stability under confinement.
            % Returns isUnstable = true if TPD < 0, indicating a phase split lowers system Gibbs energy.
            
            z = reshape(z, [], 1);
            z = max(z, 1e-20);
            z = z / sum(z);
            PV = Pext;
            
            % 1. Wilson K-factor Seeding (PV converted to MPa for correlation alignment if required)
            % Using strict thermodynamic dimensionless ratios:
            Pc_vec = obj.EOS.Fluid.Pc(:);
            Tc_vec = obj.EOS.Fluid.Tc(:);
            om_vec = obj.EOS.Fluid.omega(:);
            
            K0 = (Pc_vec ./ PV) .* exp(5.37 .* (1 + om_vec) .* (1 - Tc_vec ./ T));
            x0 = obj.projectSimplex(z ./ K0); % Liquid-like seed
            
            % 2. Baseline Vapor State Evaluation at PV
            [lnphiV, ZV, ~, ~] = obj.EOS.calculateState(PV, T, z, 1, r_cap);
            
            % 3. Initialize Liquid Pressure & Capillary Seeding
            [~, ZL0, ~, ~] = obj.EOS.calculateState(PV, T, x0, -1, r_cap);
            [rho_l0, rho_v0] = obj.calculatePhaseDensities(ZL0, ZV, PV, PV, T);
            
            Pc_seed = obj.evaluateCapillaryPressure(x0, z, rho_l0, rho_v0, r_cap);
            PL0 = max(PV - Pc_seed, 1e3); % Guard positivity
            
            % 4. Execute Projected Gradient Minimization with Lagged Pressure Coupling
            [w_star, tpd_min, PL, Pc] = obj.minimizeTPD(x0, z, T, PV, PL0, lnphiV, ZV, r_cap);
            
            % 5. Evaluate Stability Criteria & Assign Output Vectors
            isUnstable = (tpd_min < 0);
            
            if isUnstable
                epsK = 1e-14;
                w_star = obj.projectSimplex(w_star);
                K_factors = max(z ./ max(w_star, epsK), epsK); % K = y/x for dew tracing
            else
                K_factors = nan(size(z));
                Pc = max(PV - PL, 0);
            end
        end
    end
    
    methods (Access = private)
        function [w, tpd_min, PL, Pc] = minimizeTPD(obj, w0, z, T, PV, PL0, lnphiV, ZV, r_cap)
            % Core projected gradient descent engine with simplex projection
            w  = obj.projectSimplex(w0);
            PL = PL0;
            Pc = max(PV - PL0, 0);
            nc = obj.EOS.Fluid.NC;
            
            for it = 1:obj.MaxIterations
                % 1) Evaluate trial liquid state at current PL
                [lnphiL, ZL, ~, ~] = obj.EOS.calculateState(PL, T, w, -1, r_cap);
                
                % 2) Update densities and macroscopic capillary discontinuity
                [rho_l, rho_v] = obj.calculatePhaseDensities(ZL, ZV, PL, PV, T);
                
                % Trivial root protection: check similarity to feed composition
                if norm(w - z) < 1e-4
                    Pc_new = 0.0;
                else
                    Pc_new = obj.evaluateCapillaryPressure(w, z, rho_l, rho_v, r_cap);
                end
                
                if isnan(Pc_new) || Pc_new < 0, Pc_new = 0.0; end
                
                % Under-relaxed successive substitution updates for pressure gap
                if isfinite(Pc_new) && Pc_new >= 0
                    Pc = (1 - obj.RelaxPL) * Pc + obj.RelaxPL * Pc_new;
                end
                
                PL_new = max(PV - Pc, 1e3);
                if ~isfinite(PL_new), PL_new = 1e3; end
                PL = (1 - obj.RelaxPL) * PL + obj.RelaxPL * PL_new;
                
                % 3) Evaluate Central Finite-Difference Gradients
                grad = zeros(nc, 1);
                for i = 1:nc
                    e_vec = zeros(nc, 1); e_vec(i) = 1.0;
                    
                    wf = obj.projectSimplex(w + obj.Delta * e_vec);
                    wb = obj.projectSimplex(w - obj.Delta * e_vec);
                    
                    [lnphiL_f, ~, ~, ~] = obj.EOS.calculateState(PL, T, wf, -1, r_cap);
                    [lnphiL_b, ~, ~, ~] = obj.EOS.calculateState(PL, T, wb, -1, r_cap);
                    
                    tf = obj.calculateTPDValue(wf, z, lnphiL_f, lnphiV, PL, PV);
                    tb = obj.calculateTPDValue(wb, z, lnphiL_b, lnphiV, PL, PV);
                    
                    grad(i) = (tf - tb) / (2 * obj.Delta);
                end
                
                % 4) Check Gradient Norm Convergence
                if norm(grad) < obj.GradTolerance
                    break;
                end
                
                % 5) Take Projected Gradient Step
                w_new = obj.projectSimplex(w - obj.Alpha * grad);
                if norm(w_new - w, 1) < obj.StepTolerance
                    w = w_new;
                    break;
                end
                w = w_new;
            end
            
            % Final evaluation at the stationary minimiser
            [lnphiL, ~, ~, ~] = obj.EOS.calculateState(PL, T, w, -1, r_cap);
            tpd_min = obj.calculateTPDValue(w, z, lnphiL, lnphiV, PL, PV);
        end
        
        function tpd_val = calculateTPDValue(~, w, z, lnphiL, lnphiV, PL, PV)
            % Evaluates non-isobaric TPD including the mechanical pressure offset
            eps_val = 1e-12;
            chem_term = sum(w .* (log((w + eps_val) ./ z) + lnphiL - lnphiV));
            mech_term = log(max(PL, eps_val) / max(PV, eps_val));
            tpd_val   = chem_term + mech_term;
        end
        
        function Pc = evaluateCapillaryPressure(obj, x, y, rho_l, rho_v, r_cap)
            % MacLeod-Sugden Parachor IFT model coupled with Young-Laplace equation
            Pch = obj.EOS.Fluid.Parachor(:);
            
            % Molar densities in mol/cm^3 or SI conversion alignment
            % Parachor standard: IFT^(1/4) = sum( Pch_i * (x_i * rho_L - y_i * rho_V) )
            % Assuming rho_l, rho_v are in mol/m^3, convert to mol/cm^3 by dividing by 1e6
            rho_l_cgs = rho_l / 1e6;
            rho_v_cgs = rho_v / 1e6;
            
            ift_param = sum(Pch .* (x .* rho_l_cgs - y .* rho_v_cgs));
            ift_param = max(ift_param, 0.0);
            
            sigma_mN_m = ift_param^4;       % IFT in mN/m (or dynes/cm)
            sigma_N_m  = sigma_mN_m / 1000; % Convert to N/m
            
            % Young-Laplace: Pcap = (2 * sigma * cos(theta)) / r_cap
            theta_rad = deg2rad(obj.EOS.Rock.Theta);
            Pc = (2 * sigma_N_m * cos(theta_rad)) / r_cap;
        end
        
        function [rho_l, rho_v] = calculatePhaseDensities(obj, ZL, ZV, PL, PV, T)
            % Molar density evaluations [mol/m^3] via Z-factors
            R_val = 8.3144621;
            rho_l = PL / (ZL * R_val * T);
            rho_v = PV / (ZV * R_val * T);
        end
        
        function w_proj = projectSimplex(~, w)
            % Duchi et al. (2008) efficient Euclidean projection onto the probability simplex
            n = length(w);
            u = sort(w, 'descend');
            cssv = cumsum(u) - 1;
            idx = (1:n)';
            rho = find(u - cssv ./ idx > 0, 1, 'last');
            theta = cssv(rho) / rho;
            w_proj = max(w - theta, 0);
            w_proj = w_proj / sum(w_proj);
        end
        
        function parseNumericalOverrides(obj, varargin)
            % Allows overriding optimization tolerances during class construction
            for idx = 1:2:length(varargin)
                switch string(varargin{idx})
                    case "MaxIterations", obj.MaxIterations = varargin{idx+1};
                    case "GradTolerance", obj.GradTolerance = varargin{idx+1};
                    case "StepTolerance", obj.StepTolerance = varargin{idx+1};
                    case "Alpha",         obj.Alpha         = varargin{idx+1};
                    case "Delta",         obj.Delta         = varargin{idx+1};
                    case "RelaxPL",       obj.RelaxPL       = varargin{idx+1};
                end
            end
        end
    end
end