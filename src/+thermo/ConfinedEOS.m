classdef ConfinedEOS < handle
    % CONFINEDEOS Integrated Thermodynamic Confinement Equation of State Engine
    % Consolidates Peng-Robinson VLE foundations with mineralogy-driven 
    % Fluid-Wall Interactions (FWI) and Ternary Interaction Parameters (TIP).
    
    properties (SetAccess = private)
        Fluid entities.FluidProperties  % Handle to active fluid mixture data model
        Rock  entities.RockProperties   % Handle to active rock matrix mineralogy model
    end
    
    properties (Access = public)
        % Customizable Tuning Interaction Parameters (TIP framework fitting coefficients)
        k (1,1) double = 0.0            % Base confinement structural multiplier
        lambda (1,1) double = 0.0       % Acentric scaling tuning coefficient
        pT_wall (1,1) double = 1.0      % Pore-geometry rational exponent power
        kijc (:,:) double = []          % Confinement interaction parameter matrix [NC x NC]
    end
    
    properties (Constant, Access = private)
        R double = 8.3144621            % Universal Gas Constant [J/(mol*K)]
        NA double = 6.02214076e+23      % Avogadro's Constant [1/mol]
    end
    
    methods
        function obj = ConfinedEOS(fluidEntity, rockEntity, varargin)
            % Constructor binding structural models to the active thermodynamics engine
            validateattributes(fluidEntity, {'entities.FluidProperties'}, {'scalar'});
            validateattributes(rockEntity, {'entities.RockProperties'}, {'scalar'});
            obj.Fluid = fluidEntity;
            obj.Rock = rockEntity;
            
            % Enforce fallback default for confinement matrix if not explicitly assigned
            obj.kijc = zeros(obj.Fluid.NC, obj.Fluid.NC);
            
            % Parse optional tuning parameters if passed during construction
            if nargin > 2
                obj.parseOptionalTuners(varargin{:});
            end
        end
        
        function [lnphi, Z, V_shifted, confined_params] = calculateState(obj, P, T, z, phaseFlag, r_cap)
            % Main calculation entry point matching your wrapper workflow layout
            % phaseFlag: -1 for Vapor phase tracking, +1 for Liquid phase tracking
            
            z = reshape(z, [], 1); % Force strict column matrix orientation
            
            % 1. Pure Component Parameters Calculation (PR a_i, b_i and Confinement c_i)
            [a_pure, b_pure, c_pure, dH, eps_wall_mol, Tfac, size_term] = obj.calculatePureParameters(T, r_cap);
            
            % 2. Application of Overloaded Mixing Rules
            [amix_a, amix2_a]   = obj.applyMixingRule(z, a_pure, obj.Fluid.BIPMatrix);
            [amix_c, amix2_c]   = obj.applyMixingRule(z, c_pure, obj.kijc);
            
            % Consolidate effective attraction matrices
            amix_eff = max(amix_a - amix_c, 1e-28);
            amix2_eff = amix2_a - amix2_c;
            bmix = sum(z .* b_pure); % Strict linear allocation entry
            
            % 3. Extract Cubic Real Root for Compressibility (Z)
            [Z, dP] = obj.solveZFactor(P, T, amix_eff, bmix, phaseFlag);
            
            % 4. Peneloux Volume Translation
            c_shift_pure = obj.Fluid.VolumeShift' .* b_pure;
            cMix = sum(z .* c_shift_pure);
            V_ideal = (Z * obj.R * T) / P;
            V_shifted = V_ideal - cMix;
            
            % 5. Evaluate the Vectorized Analytical Fugacity Field
            [lnphi, ~] = obj.evaluateFugacityField(P, T, Z, b_pure, amix_eff, bmix, amix2_eff, dP);
            
            % 6. Populate Diagnostic Struct for Solver Verification Logs
            confined_params = struct(...
                'a_bulk_iT',  a_pure, ...
                'b',          b_pure, ...
                'c_i',        c_pure, ...
                'amix_a',     amix_a, ...
                'amix_c',     amix_c, ...
                'amix_eff',   amix_eff, ...
                'A',          (amix_eff * P) / (obj.R * T)^2, ...
                'B',          (bmix * P) / (obj.R * T), ...
                'zfactor',    Z, ...
                'rel_raw',    dH, ...
                'eps_wall_i', eps_wall_mol, ...
                'Tfac_i',     Tfac, ...
                'size_term',  size_term);
        end
    end
    
    methods (Access = private)
        function [a, b, c, dH, eps_wall_mol, Tfac, size_term] = calculatePureParameters(obj, T, r_cap)
            % Combines calcab_multicomp and calc_conf_3param logic loops
            nc = obj.Fluid.NC;
            a = zeros(nc, 1);
            b = zeros(nc, 1);
            c = zeros(nc, 1);
            
            % Calculate geometric ratio profiles
            dH = obj.Fluid.LJ_Size(:) ./ r_cap;
            mask = (dH > 0.0);
            
            % Gather normalized weight and energy foundations from entities
            w_norm = obj.Rock.WeightVector(:);
            eps_vector_mol = obj.Rock.EnergyVector(:);
            
            % Fluid energy conversion from Joules per molecule to Joules per mole
            eps_fluid_mol = obj.Fluid.LJ_Energy(:) * obj.NA;
            
            % Evaluate composite rock wall activity via geometric mean summation
            S_active = sum(w_norm .* sqrt(eps_vector_mol));
            eps_wall_mol = sqrt(eps_fluid_mol) * S_active;
            
            % Evaluate component properties array loop
            for i = 1:nc
                % Advanced Peng-Robinson Acentric factor branching thresholds
                if obj.Fluid.omega(i) > 0.49
                    m = 0.379642 + 1.48503 * obj.Fluid.omega(i) - 0.164423 * obj.Fluid.omega(i)^2 + 0.016666 * obj.Fluid.omega(i)^3;
                else
                    m = 0.37464 + 1.54226 * obj.Fluid.omega(i) - 0.26992 * obj.Fluid.omega(i)^2;
                end
                
                alpha = (1 + m * (1 - sqrt(T / obj.Fluid.Tc(i))))^2;
                a(i) = 0.457235528921 * alpha * (obj.R * obj.Fluid.Tc(i))^2 / obj.Fluid.Pc(i);
                b(i) = 0.0777960739039 * obj.R * obj.Fluid.Tc(i) / obj.Fluid.Pc(i);
            end
            
            % Confinement parameter evaluation vector loops
            k_eff = obj.k * (1 + obj.lambda * obj.Fluid.omega(:));
            Tfac = 1 - exp(-eps_wall_mol ./ (obj.R * T));
            
            term_xp = dH.^obj.pT_wall;
            size_term = term_xp ./ (1 + 4.0 * term_xp);
            
            c(mask) = k_eff(mask) .* b(mask) .* eps_wall_mol(mask) .* Tfac(mask) .* size_term(mask);
            c(~isfinite(c)) = 0;
            c = max(c, 0);
        end
        
        function [amix, amix2] = applyMixingRule(obj, z, pure_param, bip_matrix)
            % Vectorized contraction modeling reproducing calcabmix.m math
            nc = obj.Fluid.NC;
            ij_matrix = zeros(nc, nc);
            
            for i = 1:nc
                for j = 1:nc
                    ij_matrix(i,j) = sqrt(pure_param(i) * pure_param(j)) * (1 - bip_matrix(i,j));
                end
            end
            
            amix = z' * ij_matrix * z;
            amix2 = ij_matrix * z;
        end
        
        function [Z, dP] = solveZFactor(obj, P, T, amix, bmix, phaseFlag)
            % Solves for Z and returns implicit partial derivatives analytically
            A = (amix * P) / (obj.R * T)^2;
            B = (bmix * P) / (obj.R * T);
            
            c2 = B - 1;
            c1 = A - 3 * B^2 - 2 * B;
            c0 = B^3 + B^2 - A * B;
            
            zroots = roots([1, c2, c1, c0]);
            realRoots = real(zroots(abs(imag(zroots)) < 1e-10 & real(zroots) > 0));
            realRoots = sort(realRoots);
            
            if isempty(realRoots)
                error('ConfinedEOS:ExtinctionAnomaly', 'No real roots found for the cubic EOS at P=%.2f MPa.', P/1e6);
            end
            
            if phaseFlag > 0
                Z = min(realRoots);  % Liquid phase tracking preference (+1)
            else
                Z = max(realRoots);  % Vapor phase tracking preference (-1)
            end
            
            % Analytical Derivatives Evaluation (Implicit Function Theorem)
            dG_dZ = 3*Z^2 + 2*c2*Z + c1;
            dB_dP = bmix / (obj.R * T);
            dA_dP = amix / (obj.R * T)^2;
            
            dc2_dP = dB_dP;
            dc1_dP = dA_dP - 6*B*dB_dP - 2*dB_dP;
            dc0_dP = (3*B^2 + 2*B - A)*dB_dP - B*dA_dP;
            
            dG_dP = Z^2*dc2_dP + Z*dc1_dP + dc0_dP;
            
            dZdP = -dG_dP / dG_dZ;
            
            dP = {dZdP, dB_dP};
        end
        
        function [lnfugcoef, fugcoef] = evaluateFugacityField(obj, P, T, Z, b, amix, bmix, amix2, dP)
            % Vectorized analytical execution mapping calcfugcoef_multicomp.m exactly
            nc = obj.Fluid.NC;
            lnfugcoef = zeros(nc, 1);
            fugcoef = zeros(nc, 1);
            
            c1 = 1 + sqrt(2);
            c2 = 1 - sqrt(2);
            c0 = c2 - c1;
            
            A = (amix * P) / (obj.R * T)^2;
            B = (bmix * P) / (obj.R * T);
            
            if Z < B
                error('ConfinedEOS:InvalidZFactor', 'Z-factor (%.6f) must be larger than Bmix (%.6f) to prevent logarithmic singularities.', Z, B);
            end
            
            f = log((Z + c2*B)/(Z + c1*B));
            
            for i = 1:nc
                term1 = (b(i)/bmix)*(Z - 1) - log(Z - B);
                term2 = (A / (c0 * B)) * (2 * amix2(i) / amix - b(i) / bmix) * f;
                
                lnfugcoef(i) = term1 - term2;
                fugcoef(i) = exp(lnfugcoef(i));
            end
        end
        
        function parseOptionalTuners(obj, varargin)
            % Processes name-value input overrides dynamically
            for idx = 1:2:length(varargin)
                switch string(varargin{idx})
                    case "k",       obj.k = varargin{idx+1};
                    case "lambda",  obj.lambda = varargin{idx+1};
                    case "pT_wall", obj.pT_wall = varargin{idx+1};
                    case "kijc",    obj.kijc = varargin{idx+1};
                end
            end
        end
    end
end