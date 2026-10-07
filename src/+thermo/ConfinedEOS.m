classdef ConfinedEOS < handle
    % CONFINEDEOS Peng-Robinson EOS with a mineralogy-driven fluid-wall
    % interaction (FWI) attraction correction:
    %
    %   a_eff,ij = sqrt(a_i a_j)(1 - k_ij) - sqrt(c_i c_j)(1 - kc_ij)
    %   c_i      = k (1 + lambda*omega_i) * b_i * W_i * s(sigma_i / r)
    %   s(x)     = x^p / (1 + 4 x^p)
    %   W_i      = eps_fw,i * (1 - exp(-eps_fw,i / RT))           ("meanfield")
    %            = sum_m phi_m eps_im (1 - exp(-eps_im / RT))     ("patch")
    %   eps_fw,i = sqrt(eps_f,i) * sum_m phi_m sqrt(eps_m)
    %   eps_im   = sqrt(eps_f,i * eps_m)
    %
    % phi_m are per-mineral mass or volume fractions from RockProperties
    % (WallWeighting). c_i is composition-independent, so the quadratic
    % mixing rule keeps the standard PR fugacity expression exact.
    %
    % Units: P [Pa], T [K], r_cap [m], energies [J/mol], a [Pa m^6/mol^2],
    % b [m^3/mol].

    properties (SetAccess = private)
        Fluid entities.FluidProperties
        Rock  entities.RockProperties
    end

    properties (Access = public)
        % FWI tuning parameters (fitted on literature Tc-shift data)
        k       (1,1) double {mustBeReal, mustBeFinite} = 64.4376
        lambda  (1,1) double {mustBeReal, mustBeFinite} = 10.0338
        pT_wall (1,1) double {mustBeReal, mustBeFinite} = 1.3745
        kijc    (:,:) double {mustBeReal, mustBeFinite} = []   % [NC x NC], symmetric

        % Homotopy knob: c_i -> ConfinementScale * c_i. 1 = full confinement,
        % 0 = bulk. Used to continue a dew point from the bulk solution.
        ConfinementScale (1,1) double {mustBeInRange(ConfinementScale, 0, 1)} = 1

        % Model options
        WallWeighting (1,1) string {mustBeMember(WallWeighting, ["mass","volume"])} = "mass"
        WallMixing    (1,1) string {mustBeMember(WallMixing, ["meanfield","patch"])} = "meanfield"
        ShiftFugacity (1,1) logical = true   % apply Peneloux shift to ln(phi)


    end

    properties (Constant, Access = private)
        R  = 8.3144621          % [J/(mol K)] (same value as RockProperties)
        NA = 6.02214076e+23     % [1/mol]
        OmegaA = 0.457235528921
        OmegaB = 0.0777960739039
        d1 = 1 + sqrt(2)
        d2 = 1 - sqrt(2)
        RootImagTol = 1e-10
    end

    methods
        function obj = ConfinedEOS(fluidEntity, rockEntity, opts)
            arguments
                fluidEntity (1,1) entities.FluidProperties
                rockEntity  (1,1) entities.RockProperties
                opts.k       (1,1) double
                opts.lambda  (1,1) double
                opts.pT_wall (1,1) double
                opts.kijc    (:,:) double
                opts.WallWeighting (1,1) string
                opts.WallMixing    (1,1) string
                opts.ShiftFugacity (1,1) logical
                opts.ConfinementScale (1,1) double
            end
            if fluidEntity.NC == 0
                error('ConfinedEOS:EmptyFluid', 'FluidProperties contains no components.');
            end
            if rockEntity.NumMinerals == 0
                error('ConfinedEOS:EmptyRock', 'RockProperties contains no minerals.');
            end

            obj.Fluid = fluidEntity;
            obj.Rock  = rockEntity;
            obj.kijc  = zeros(fluidEntity.NC);

            % Unknown names are rejected by the arguments block (no silent typos)
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function set.kijc(obj, K)
            if ~isempty(obj.Fluid) && ~isempty(K)  %#ok<MCSUP>
                n = obj.Fluid.NC;                  %#ok<MCSUP>
                if ~isequal(size(K), [n n])
                    error('ConfinedEOS:KijcSize', 'kijc must be %dx%d.', n, n);
                end
                if max(abs(K - K.'), [], 'all') > 1e-12
                    error('ConfinedEOS:KijcAsymmetric', 'kijc must be symmetric.');
                end
            end
            obj.kijc = K;
        end

        function [lnphi, Z, V_shifted, confined_params] = calculateState(obj, P, T, z, phaseFlag, r_cap)
            % P [Pa], T [K], z mole fractions, phaseFlag: +1 liquid root,
            % -1 vapor root, 0 minimum-Gibbs root; r_cap [m] (Inf = bulk).
            arguments
                obj
                P (1,1) double {mustBePositive, mustBeFinite}
                T (1,1) double {mustBePositive, mustBeFinite}
                z (:,1) double {mustBeNonnegative}
                phaseFlag (1,1) double {mustBeMember(phaseFlag, [-1 0 1])}
                r_cap (1,1) double {mustBePositive}
            end
            if numel(z) ~= obj.Fluid.NC
                error('ConfinedEOS:CompositionSize', 'z has %d entries; NC = %d.', numel(z), obj.Fluid.NC);
            end
            if abs(sum(z) - 1) > 1e-8
                error('ConfinedEOS:CompositionSum', 'Mole fractions sum to %.10f.', sum(z));
            end
            if isfinite(r_cap) && r_cap < 1e-10
                error('ConfinedEOS:PoreRadiusUnits', ...
                    'r_cap = %.3g looks non-SI; pass the pore radius in metres.', r_cap);
            end

            RT = obj.R * T;

            % 1. Pure-component parameters
            [a, b, c, dH, eps_fw, Tfac, size_term] = obj.calculatePureParameters(T, r_cap);

            % 2. Mixing rules (bulk attraction minus confinement term)
            [amix_a, Sa] = obj.quadraticMix(z, a, obj.Fluid.BIPMatrix);
            [amix_c, Sc] = obj.confinementMixing(z, c, eps_fw, RT);
            amix = amix_a - amix_c;
            Si   = Sa - Sc;          % (1/2n) d(n^2 a_eff)/dn_i ; Sc = c_hat_i (= (Cx)_i if beta = 0)
            bmix = z.' * b;

            if amix <= 0
                error('ConfinedEOS:NonPositiveAttraction', ...
                    ['Effective attraction a_eff = %.4g <= 0 (confinement term exceeds ' ...
                     'bulk attraction) at r = %.3g m.'], amix, r_cap);
            end

            A = amix * P / RT^2;
            B = bmix * P / RT;

            % 3. Compressibility factor
            [Z, nRoots] = obj.solveZFactor(A, B, phaseFlag);

            % 4. Fugacity coefficients (PR, quadratic mixing rule)
            lnphi = (b / bmix) * (Z - 1) - log(Z - B) ...
                  - A / (2*sqrt(2)*B) * (2 * Si / amix - b / bmix) ...
                    * log((Z + obj.d1*B) / (Z + obj.d2*B));

            % 5. Peneloux volume translation
            c_shift = obj.Fluid.VolumeShift(:) .* b;
            V_shifted = Z * RT / P - z.' * c_shift;
            if obj.ShiftFugacity
                lnphi = lnphi - c_shift * P / RT;
            end

            confined_params = struct( ...
                'a_bulk_iT',  a, ...
                'b',          b, ...
                'c_i',        c, ...
                'amix_a',     amix_a, ...
                'amix_c',     amix_c, ...
                'amix_eff',   amix, ...
                'A',          A, ...
                'B',          B, ...
                'zfactor',    Z, ...
                'nRoots',     nRoots, ...
                'rel_raw',    dH, ...
                'eps_wall_i', eps_fw, ...
                'Tfac_i',     Tfac, ...
                'size_term',  size_term, ...
                'c_shift_i',  c_shift, ...
                'c_hat_i',    Sc, ...
                'WallWeighting', obj.WallWeighting, ...
                'WallMixing',    obj.WallMixing, ...
                'ConfinementScale', obj.ConfinementScale);
        end

        function [cm, c_hat] = confinementMix(obj, x, T, r_cap)
            % Mixture confinement term and its partial-molar derivative at
            % composition x (public access for diagnostics).
            [~, ~, c] = obj.calculatePureParameters(T, r_cap);
            [cm, c_hat] = obj.confinementMixing(x(:), c, [], []);
        end

        function lockWall(obj, rockEntity)
            % Replace the wall description (e.g. a single-mineral calibration
            % wall) while keeping the fluid, TIPs and tuning parameters.
            arguments
                obj
                rockEntity (1,1) entities.RockProperties
            end
            if rockEntity.NumMinerals == 0
                error('ConfinedEOS:EmptyRock', 'RockProperties contains no minerals.');
            end
            obj.Rock = rockEntity;
        end

        function [a, b, c] = pureParameters(obj, T, r_cap)
            % Pure-component a_i(T), b_i and confinement c_i(T, r) [column vectors].
            % Public read-only access for calibration and diagnostics.
            arguments
                obj
                T     (1,1) double {mustBePositive}
                r_cap (1,1) double {mustBePositive}
            end
            [a, b, c] = obj.calculatePureParameters(T, r_cap);
        end

        function [Tcp, Pcp] = pureCriticalPoint(obj, i, r_cap)
            % Confined critical point of pure component i at pore radius r_cap [m].
            % With b fixed and a_eff(T) = a_i(T) - c_i(T) independent of volume,
            % dP/dV = d2P/dV2 = 0 reduce to the PR conditions evaluated with a_eff:
            %   a_eff(Tcp) / (b R Tcp) = OmegaA / OmegaB,   Pcp = OmegaB R Tcp / b
            % Consequence: the model predicts dPc/Pc = dTc/Tc (b is unconfined).
            arguments
                obj
                i     (1,1) double {mustBeInteger, mustBePositive}
                r_cap (1,1) double {mustBePositive}
            end
            Tc = obj.Fluid.Tc(i);
            [~, b, c] = obj.calculatePureParameters(Tc, r_cap);
            if c(i) == 0
                Tcp = Tc;
            else
                f   = @(T) obj.criticalResidual(T, i, r_cap);
                Thi = Tc;                              % f(Tc) = -c/(bRTc) < 0
                Tlo = 0.9 * Tc;
                while f(Tlo) <= 0
                    Thi = Tlo;
                    Tlo = 0.9 * Tlo;
                    if Tlo < 0.05 * Tc
                        error('ConfinedEOS:NoConfinedCriticalPoint', ...
                            'No confined critical point above 0.05 Tc for component %d.', i);
                    end
                end
                Tcp = fzero(f, [Tlo, Thi], optimset('TolX', 1e-13 * Tc));   % tight: Stage-1 gradients are finite differences of this root
            end
            Pcp = obj.OmegaB * obj.R * Tcp / b(i);
        end

        function Pcap = capillaryPressure(obj, x, y, V_L, V_V, r_cap)
            % Young-Laplace with MacLeod-Sugden IFT. Single source for all solvers.
            % x, y : liquid / vapor mole fractions; V_L, V_V : Peneloux-translated
            % molar volumes [m^3/mol] (the parachor correlation is calibrated on
            % real densities); r_cap [m]. Returns Pcap [Pa] (>= 0 for theta < 90 deg).
            if ~isfinite(r_cap)
                Pcap = 0.0;
                return;
            end
            rhoL = 1e-6 / V_L;                      % [mol/cm^3]
            rhoV = 1e-6 / V_V;
            s14  = obj.Fluid.Parachor(:).' * (x(:) * rhoL - y(:) * rhoV);
            sigma = max(s14, 0)^4 * 1e-3;           % [mN/m] -> [N/m]
            Pcap = 2 * sigma * cosd(obj.Rock.Theta) / r_cap;
        end

        function K = wilsonK(obj, P, T)
            % Wilson K-factor estimate (single definition for all solvers).
            K = (obj.Fluid.Pc(:) / P) .* exp(5.373 * (1 + obj.Fluid.omega(:)) .* (1 - obj.Fluid.Tc(:) / T));
        end

        function [eps_fw, W, Tfac] = wallEnergy(obj, T)
            % Per-component fluid-wall energy [J/mol], W_i = <eps (1 - e^{-eps/RT})>,
            % and effective Tfac = W ./ eps_fw.
            RT = obj.R * T;
            switch obj.WallWeighting
                case "mass",   phi = obj.Rock.MassFraction;
                case "volume", phi = obj.Rock.VolumeFraction;
            end
            present = phi > 0;
            phi = phi(present);
            Em  = obj.Rock.Energy_Jmol(present);            % 1 x M
            eps_f = obj.Fluid.LJ_Energy(:) * obj.NA;        % nc x 1 [J/mol]

            switch obj.WallMixing
                case "meanfield"
                    eps_fw = sqrt(eps_f) * (phi * sqrt(Em).');
                    Tfac   = 1 - exp(-eps_fw / RT);
                    W      = eps_fw .* Tfac;
                case "patch"
                    eps_im = sqrt(eps_f * Em);              % nc x M
                    W      = (eps_im .* (1 - exp(-eps_im / RT))) * phi.';
                    eps_fw = eps_im * phi.';
                    Tfac   = W ./ eps_fw;
            end
        end
    end

    methods (Access = private)
        function [a, b, c, dH, eps_fw, Tfac, size_term] = calculatePureParameters(obj, T, r_cap)
            Tc = obj.Fluid.Tc(:);
            Pc = obj.Fluid.Pc(:);
            w  = obj.Fluid.omega(:);

            % PR78 alpha function
            m = 0.37464 + 1.54226*w - 0.26992*w.^2;
            hi = w > 0.49;
            m(hi) = 0.379642 + 1.48503*w(hi) - 0.164423*w(hi).^2 + 0.016666*w(hi).^3;
            alpha = (1 + m .* (1 - sqrt(T ./ Tc))).^2;

            a = obj.OmegaA * alpha .* (obj.R * Tc).^2 ./ Pc;
            b = obj.OmegaB * obj.R * Tc ./ Pc;

            % Confinement term (r_cap = Inf -> dH = 0 -> c = 0, bulk limit)
            dH = obj.Fluid.LJ_Size(:) / r_cap;
            xp = dH .^ obj.pT_wall;
            size_term = xp ./ (1 + 4*xp);

            [eps_fw, W, Tfac] = obj.wallEnergy(T);
            k_eff = obj.k * (1 + obj.lambda * w);

            c = obj.ConfinementScale * (k_eff .* b .* W .* size_term);
            if any(~isfinite(c)) || any(c < 0)
                error('ConfinedEOS:InvalidConfinementTerm', ...
                    'Non-finite or negative c_i (check k, lambda, pT_wall and r_cap).');
            end
        end

        function g = criticalResidual(obj, T, i, r_cap)
            [a, b, c] = obj.calculatePureParameters(T, r_cap);
            g = (a(i) - c(i)) / (b(i) * obj.R * T) - obj.OmegaA / obj.OmegaB;
        end

        function [cm, c_hat] = confinementMixing(obj, x, c, ~, ~)
            % Mole-fraction quadratic rule, identical in form to a_mix:
            %   c_m = sum_i sum_j x_i x_j sqrt(c_i c_j) (1 - TIP_ij)
            %   c_hat_i = (1/2n) d(n^2 c_m)/dn_i = (C x)_i
            % (The surface-fraction variant was tested and rejected.)
            [cm, c_hat] = obj.quadraticMix(x, c, obj.kijc);
        end

        function [Z, nRoots] = solveZFactor(obj, A, B, phaseFlag)
            zr = roots([1, B - 1, A - 3*B^2 - 2*B, B^3 + B^2 - A*B]);
            isReal = abs(imag(zr)) <= obj.RootImagTol * max(1, abs(real(zr)));
            zr = sort(real(zr(isReal)));
            zr = zr(zr > B);                               % physical branch only
            nRoots = numel(zr);
            if nRoots == 0
                error('ConfinedEOS:NoPhysicalRoot', 'No real root with Z > B (A = %.4g, B = %.4g).', A, B);
            end

            switch phaseFlag
                case  1, Z = zr(1);
                case -1, Z = zr(end);
                case  0
                    % Dimensionless residual Gibbs energy sum_i z_i ln(phi_i)
                    g = zr - 1 - log(zr - B) ...
                        - A/(2*sqrt(2)*B) * log((zr + obj.d1*B) ./ (zr + obj.d2*B));
                    [~, iMin] = min(g);
                    Z = zr(iMin);
            end
        end
    end

    methods (Static, Access = private)
        function [amix, Si] = quadraticMix(z, p, K)
            Aij  = sqrt(p * p.') .* (1 - K);
            Si   = Aij * z;
            amix = z.' * Si;
        end
    end
end
