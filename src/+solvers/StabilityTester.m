classdef StabilityTester < handle
    % STABILITYTESTER Non-isobaric & Isobaric Tangent Plane Distance (TPD) stability engine.
    % Evaluates phase stability under nanoporous confinement by coupling projected
    % gradient descent on the probability simplex with optional lagged capillary pressure updates.
    % Supports both Combined (FWI + Pcap) and FWI-Only (Pcap = 0) physics modes.

    properties (SetAccess = private)
        EOS thermo.ConfinedEOS  % Handle to active thermodynamic operator
        MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 500
        GradTolerance (1,1) double {mustBePositive, mustBeReal} = 1e-6
        StepTolerance (1,1) double {mustBePositive, mustBeReal} = 1e-6
        Alpha (1,1) double {mustBePositive, mustBeReal} = 1e-2       % Projected gradient step size
        Delta (1,1) double {mustBePositive, mustBeReal} = 1e-6       % Finite-difference perturbation step
        RelaxPL (1,1) double {mustBeReal} = 0.5                      % Under-relaxation factor for pressure updates
    end
    
    methods
        function obj = StabilityTester(eosEngine, varargin)
            % Constructor attaching physical models and numerical boundaries
            validateattributes(eosEngine, {'thermo.ConfinedEOS'}, {'scalar'});
            obj.EOS = eosEngine;
            
            % Process input overrides if passed during construction
            if nargin > 1
                obj.parseNumericalOverrides(varargin{:});
            end
        end
        
        function [isUnstable, w_star, tpd_min, K_factors, PV, PL, Pc] = evaluateDewStability(obj, Pext, T, z, r_cap, varargin)
            % Main entry point evaluating vapor-phase (dew) stability under confinement.
            % Returns isUnstable = true if TPD < 0, indicating a phase split lowers system Gibbs energy.
            % Optional Key-Value input:
            %   'CapMode': 'Pc' (Combined Model) or 'PL' / 'NoCap' (FWI-Only Model)
            
            z = reshape(z, [], 1);
            z = max(z, 1e-20);
            z = z / sum(z);
            PV = Pext;
            
            % Parse Mode Controls
            capMode = 'Pc'; % Default to Combined Model
            for idx = 1:2:length(varargin)
                if strcmpi(string(varargin{idx}), "CapMode")
                    capMode = varargin{idx+1};
                end
            end
            
            isCapActive = ~strcmpi(capMode, 'PL') && ~strcmpi(capMode, 'NoCap') && ~isinf(r_cap);

            % 1. Wilson K-factor Seeding
            Pc_vec = obj.EOS.Fluid.Pc(:);
            Tc_vec = obj.EOS.Fluid.Tc(:);
            om_vec = obj.EOS.Fluid.omega(:);

            K0 = (Pc_vec ./ PV) .* exp(5.37 .* (1 + om_vec) .* (1 - Tc_vec ./ T));
            x0 = obj.projectSimplex(z ./ K0); % Liquid-like seed

            % 2. Baseline Vapor State Evaluation at PV
            [lnphiV, ZV, ~, ~] = obj.EOS.calculateState(PV, T, z, -1, r_cap);

            % 3. Initialize Liquid Pressure & Capillary Seeding (+1 = Liquid)
            if isCapActive
                [~, ZL0, ~, ~] = obj.EOS.calculateState(PV, T, x0, 1, r_cap);
                [rho_l0, rho_v0] = obj.calculatePhaseDensities(ZL0, ZV, PV, PV, T);
                Pc_seed = obj.evaluateCapillaryPressure(x0, z, rho_l0, rho_v0, r_cap);
                PL0 = max(PV - Pc_seed, 1e3);
            else
                PL0 = PV;
                Pc_seed = 0.0;
            end
            
            % 4. Execute Projected Gradient Minimization
            [w_star, tpd_min, PL, Pc] = obj.minimizeTPD(x0, z, T, PV, PL0, lnphiV, ZV, r_cap, isCapActive);
            
            % 5. Evaluate Stability Criteria & Assign Output Vectors
            isUnstable = (tpd_min < -1e-7);
            
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
        function [w, tpd_min, PL, Pc] = minimizeTPD(obj, w0, z, T, PV, PL0, lnphiV, ZV, r_cap, isCapActive)
            % Projected gradient descent engine supporting both isobaric (FWI-Only) 
            % and non-isobaric (Combined) TPD minimization.

            w  = obj.projectSimplex(w0);
            PL = PL0;
            Pc = max(PV - PL0, 0);
            nc = obj.EOS.Fluid.NC;
            tpd_old = 1e10;

            for it = 1:obj.MaxIterations
                % 1) Evaluate trial liquid state at current PL
                [lnphiL, ZL, ~, ~] = obj.EOS.calculateState(PL, T, w, 1, r_cap);
                
                % 2) Calculate current TPD Value
                tpd_val = obj.calculateTPDValue(w, z, lnphiL, lnphiV, PL, PV, isCapActive);
                
                % Check convergence on functional change (from FWI-Only legacy logic)
                if abs(tpd_val - tpd_old) < obj.StepTolerance && it > 1
                    break;
                end
                tpd_old = tpd_val;
                
                % 3) Update capillary pressure and PL if capillary mode is active
                if isCapActive
                    [rho_l, rho_v] = obj.calculatePhaseDensities(ZL, ZV, PL, PV, T);
                    
                    if norm(w - z) < 1e-4
                        Pc_new = 0.0;
                    else
                        Pc_new = obj.evaluateCapillaryPressure(w, z, rho_l, rho_v, r_cap);
                    end
                    
                    if isnan(Pc_new) || Pc_new < 0, Pc_new = 0.0; end
                    Pc = (1 - obj.RelaxPL) * Pc + obj.RelaxPL * Pc_new;
                    
                    PL_new = max(PV - Pc, 1e3);
                    PL = (1 - obj.RelaxPL) * PL + obj.RelaxPL * PL_new;
                else
                    PL = PV;
                    Pc = 0.0;
                end
                
                % 4) Evaluate Finite-Difference Gradients
                grad = zeros(nc, 1);
                for i = 1:nc
                    e_vec = zeros(nc, 1); e_vec(i) = 1.0;
                    wf = obj.projectSimplex(w + obj.Delta * e_vec);
                    wb = obj.projectSimplex(w - obj.Delta * e_vec);
                    
                    [lnphiL_f, ~, ~, ~] = obj.EOS.calculateState(PL, T, wf, 1, r_cap);
                    [lnphiL_b, ~, ~, ~] = obj.EOS.calculateState(PL, T, wb, 1, r_cap);
                    
                    tf = obj.calculateTPDValue(wf, z, lnphiL_f, lnphiV, PL, PV, isCapActive);
                    tb = obj.calculateTPDValue(wb, z, lnphiL_b, lnphiV, PL, PV, isCapActive);
                    grad(i) = (tf - tb) / (2 * obj.Delta);
                end
                
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
            
            % Final TPD evaluation at stationary point
            [lnphiL, ~, ~, ~] = obj.EOS.calculateState(PL, T, w, 1, r_cap);
            tpd_min = obj.calculateTPDValue(w, z, lnphiL, lnphiV, PL, PV, isCapActive);
        end
        
        function tpd_val = calculateTPDValue(~, w, z, lnphiL, lnphiV, PL, PV, isCapActive)
            eps_val = 1e-12;
            chem_term = sum(w .* (log((w + eps_val) ./ z) + lnphiL - lnphiV));
            
            if isCapActive
                mech_term = log(max(PL, eps_val) / max(PV, eps_val));
            else
                mech_term = 0.0; % FWI-Only Mode has no mechanical pressure jump
            end
            
            tpd_val = chem_term + mech_term;
        end
        
        function Pc = evaluateCapillaryPressure(obj, x, y, rho_l, rho_v, r_cap)
            % MacLeod-Sugden Parachor IFT model coupled with Young-Laplace equation
            Pch = obj.EOS.Fluid.Parachor(:);
            
            % Convert molar densities from mol/m^3 to gmol/cm^3
            rho_l_cgs = rho_l / 1e6;
            rho_v_cgs = rho_v / 1e6;
            
            ift_param = sum(Pch .* (x .* rho_l_cgs - y .* rho_v_cgs));
            sigma_mN_m = max(ift_param, 0.0)^4;
            sigma_N_m  = sigma_mN_m / 1000; % Convert to N/m
            
            theta_rad = deg2rad(obj.EOS.Rock.Theta);
            Pc = (2 * sigma_N_m * cos(theta_rad)) / r_cap;
        end
        
        function [rho_l, rho_v] = calculatePhaseDensities(obj, ZL, ZV, PL, PV, T)
            % Molar density evaluations [mol/m^3] via Z-factors
            R_val = 8.3144621;
            rho_l = PL / (ZL * R_val * T);
            rho_v = PV / (ZV * R_val * T);
        end
        
        function w_proj = projectFeasibleDomain(~, w)
            % Redundant proxy method mapping interface standard extensions if needed
            w_proj = w;
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
                switch lower(string(varargin{idx}))
                    case "maxiterations", obj.MaxIterations = varargin{idx+1};
                    case "gradtolerance", obj.GradTolerance = varargin{idx+1};
                    case "steptolerance", obj.StepTolerance = varargin{idx+1};
                    case "alpha",         obj.Alpha         = varargin{idx+1};
                    case "delta",         obj.Delta         = varargin{idx+1};
                    case "relaxpl",       obj.RelaxPL       = varargin{idx+1};
                end
            end
        end
    end
end