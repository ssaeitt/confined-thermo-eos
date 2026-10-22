classdef EquivalentPR
    % EQUIVALENTPR Exact standard-Peng-Robinson representation of the confined
    % EOS at a fixed temperature, pore radius and rock (FWI-only).
    %
    % ConfinedEOS uses
    %   a_eff,ij = sqrt(a_i a_j)(1 - BIP_ij) - sqrt(c_i c_j)(1 - TIP_ij)
    % with composition-independent c_i(T, r) and b_i unchanged. Defining
    %   a'_i      = a_i - c_i                        -> Omega_a,i' = Omega_a (1 - c_i/a_i)
    %   BIP'_ij   = 1 - a_eff,ij / sqrt(a'_i a'_j)
    % gives sqrt(a'_i a'_j)(1 - BIP'_ij) = a_eff,ij term by term, so a standard
    % PR with per-component Omega_a' (Omega_b, Tc, Pc, omega, volume shift
    % unchanged) reproduces Z, ln phi and every phase-equilibrium result of
    % the confined model EXACTLY, at this (T, r_cap, rock). No regression.
    %
    % Validity:
    %   - Isothermal only: c_i(T) and the PR alpha differ in T, so Omega_a'
    %     and BIP' are valid at the T they were built for (T-specific deck).
    %   - FWI-only: the Combined (Pcap) mode needs a two-pressure flash and has
    %     no standard-PR equivalent.
    %   - Requires a'_i > 0 for every component (checked).
    %   - Tc, Pc, Vc are unchanged, so an LBC viscosity computed by the
    %     simulator remains the BULK viscosity correlation (stated assumption).
    %   - The simulator must apply Peneloux shifts with dimensionless s_i = c_i/b_i
    %     (FluidProperties.VolumeShift convention).
    %
    % Usage:
    %   eq = thermo.EquivalentPR.fromEOS(eos, 293.15, 11.25e-9, "B1-11nm");
    %   err = eq.verify(eos, Plist, zlist);         % max |d ln phi|, |dZ|
    %   eq.writeTable("results/simulator_eos/B1.csv");
    %   thermo.EquivalentPR.writeE300(eqArray, "results/simulator_eos/props_E300.inc");
    %   thermo.EquivalentPR.writeGEM(eqArray,  "results/simulator_eos/eos_GEM.inc");

    properties (Constant)
        OmegaA0 = 0.457235528921
        OmegaB0 = 0.0777960739039
        R = 8.3144621
    end

    properties (SetAccess = private)
        Label          (1,1) string = ""
        T              (1,1) double = NaN      % [K]
        r_cap          (1,1) double = NaN      % [m]
        RockName       (1,1) string = ""
        ComponentNames (1,:) string = string.empty(1,0)
        Tc    (1,:) double                     % [K]
        Pc    (1,:) double                     % [Pa]
        omega (1,:) double
        MW    (1,:) double
        Vc    (1,:) double                     % as stored in FluidProperties
        VolumeShift (1,:) double               % dimensionless s_i
        Parachor    (1,:) double
        OmegaA (1,:) double                    % per component, confined
        OmegaB (1,:) double                    % per component (unchanged)
        BIP    (:,:) double                    % effective, symmetric, zero diagonal
        CoverA (1,:) double                    % c_i / a_i diagnostic
    end

    methods (Static)
        function obj = fromEOS(eos, T, r_cap, label)
            arguments
                eos   (1,1) thermo.ConfinedEOS
                T     (1,1) double {mustBePositive}
                r_cap (1,1) double {mustBePositive}
                label (1,1) string = ""
            end
            fl = eos.Fluid;
            [a, ~, c] = eos.pureParameters(T, r_cap);
            ae = a - c;
            if any(ae <= 0)
                error('EquivalentPR:NonPositiveAttraction', ...
                    'a_i - c_i <= 0 for: %s (no PR equivalent at r = %.3g m).', ...
                    strjoin(fl.ComponentNames(ae <= 0), ', '), r_cap);
            end
            K = eos.kijc;
            if isempty(K), K = zeros(fl.NC); end
            aij = sqrt(a * a.') .* (1 - fl.BIPMatrix) - sqrt(c * c.') .* (1 - K);
            BIPe = 1 - aij ./ sqrt(ae * ae.');
            BIPe(1:fl.NC+1:end) = 0;
            BIPe = (BIPe + BIPe.') / 2;

            obj = thermo.EquivalentPR();
            obj.Label = label;
            obj.T = T;
            obj.r_cap = r_cap;
            obj.RockName = eos.Rock.RockName;
            obj.ComponentNames = fl.ComponentNames;
            obj.Tc = fl.Tc;  obj.Pc = fl.Pc;  obj.omega = fl.omega;  obj.MW = fl.MW;
            obj.Vc = fl.Vc;  obj.VolumeShift = fl.VolumeShift;  obj.Parachor = fl.Parachor;
            obj.OmegaA = thermo.EquivalentPR.OmegaA0 * (ae ./ a).';
            obj.OmegaB = thermo.EquivalentPR.OmegaB0 * ones(1, fl.NC);
            obj.BIP = BIPe;
            obj.CoverA = (c ./ a).';
        end

        function writeE300(eqs, filePath)
            % ECLIPSE 300 PROPS include, METRIC units (K, barsa), one record per
            % EOS region in the order of eqs (EOSNUM = 1..numel(eqs)). Requires
            % TABDIMS item 9 (NTEOS) >= numel(eqs). Verify keyword syntax and
            % units against the manual of your simulator version before use.
            arguments
                eqs (1,:) thermo.EquivalentPR
                filePath (1,1) string
            end
            thermo.EquivalentPR.checkCompatible(eqs);
            fid = thermo.EquivalentPR.openFile(filePath);
            c = onCleanup(@() fclose(fid));
            e1 = eqs(1);
            fprintf(fid, '-- Confined EOS, exact standard-PR equivalent (FWI-only), T = %.2f K\n', e1.T);
            for k = 1:numel(eqs)
                fprintf(fid, '-- EOSNUM %d : %s | rock %s | r_p = %.3f nm\n', k, eqs(k).Label, ...
                        eqs(k).RockName, eqs(k).r_cap * 1e9);
            end
            fprintf(fid, '-- Generated %s. Valid ONLY isothermally at this T.\n', ...
                    string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
            fprintf(fid, '-- Set RTEMP (or the temperature table) to %.2f degC in the deck.\n\n', e1.T - 273.15);
            fprintf(fid, 'EOS\n');
            for k = 1:numel(eqs), fprintf(fid, 'PR /\n'); end
            fprintf(fid, '\nPRCORR\n\nCNAMES\n');
            fprintf(fid, '%s\n', e1.ComponentNames);
            fprintf(fid, '/\n');
            thermo.EquivalentPR.e300Vec(fid, 'TCRIT',  eqs, @(e) e.Tc);
            thermo.EquivalentPR.e300Vec(fid, 'PCRIT',  eqs, @(e) e.Pc / 1e5);
            thermo.EquivalentPR.e300Vec(fid, 'ACF',    eqs, @(e) e.omega);
            thermo.EquivalentPR.e300Vec(fid, 'MW',     eqs, @(e) e.MW);
            thermo.EquivalentPR.e300Vec(fid, 'OMEGAA', eqs, @(e) e.OmegaA);
            thermo.EquivalentPR.e300Vec(fid, 'OMEGAB', eqs, @(e) e.OmegaB);
            thermo.EquivalentPR.e300Vec(fid, 'SSHIFT', eqs, @(e) e.VolumeShift);
            fprintf(fid, '\nBIC\n');
            for k = 1:numel(eqs)
                B = eqs(k).BIP;
                for i = 2:size(B, 1)
                    fprintf(fid, ' %.10f', B(i, 1:i-1));
                    fprintf(fid, '\n');
                end
                fprintf(fid, '/\n');
            end
        end

        function writeGEM(eqs, filePath)
            % CMG GEM fluid-model block, one EOS set per element of eqs
            % (assign with the grid EOS-type keyword). GEM expects Pc in atm.
            % Verify keyword names/units against the GEM manual of your
            % version before use; this is a template, not a validated deck.
            arguments
                eqs (1,:) thermo.EquivalentPR
                filePath (1,1) string
            end
            thermo.EquivalentPR.checkCompatible(eqs);
            fid = thermo.EquivalentPR.openFile(filePath);
            c = onCleanup(@() fclose(fid));
            e1 = eqs(1);
            nc = numel(e1.ComponentNames);
            fprintf(fid, '** Confined EOS, exact standard-PR equivalent (FWI-only), T = %.2f K\n', e1.T);
            fprintf(fid, '** Valid ONLY isothermally at this T. Pc in atm, Tc in K.\n');
            fprintf(fid, '*MODEL *PR\n*NC %d %d\n', nc, nc);
            fprintf(fid, '*COMPNAME');
            fprintf(fid, ' ''%s''', e1.ComponentNames);
            fprintf(fid, '\n');
            for k = 1:numel(eqs)
                e = eqs(k);
                fprintf(fid, '\n** ---- EOS set %d : %s | rock %s | r_p = %.3f nm\n', k, e.Label, ...
                        e.RockName, e.r_cap * 1e9);
                if numel(eqs) > 1, fprintf(fid, '*EOSSET %d\n', k); end
                thermo.EquivalentPR.gemVec(fid, '*TCRIT', e.Tc);
                thermo.EquivalentPR.gemVec(fid, '*PCRIT', e.Pc / 101325);
                thermo.EquivalentPR.gemVec(fid, '*AC',    e.omega);
                thermo.EquivalentPR.gemVec(fid, '*MW',    e.MW);
                thermo.EquivalentPR.gemVec(fid, '*OMEGA', e.OmegaA);
                thermo.EquivalentPR.gemVec(fid, '*OMEGB', e.OmegaB);
                thermo.EquivalentPR.gemVec(fid, '*VSHIFT', e.VolumeShift);
                thermo.EquivalentPR.gemVec(fid, '*PCHOR', e.Parachor);
                fprintf(fid, '*BIN\n');
                for i = 2:nc
                    fprintf(fid, ' %.10f', e.BIP(i, 1:i-1));
                    fprintf(fid, '\n');
                end
            end
            fprintf(fid, '\n*TRES %.2f   ** deg C\n', e1.T - 273.15);
        end
    end

    methods
        function [errLnPhi, errZ] = verify(obj, eos, Plist, zlist)
            % Max |ln phi_conf - ln phi_PR'| and |Z_conf - Z_PR'| over pressures
            % [Pa] and compositions (columns of zlist), for liquid, vapour and
            % minimum-Gibbs roots. Must be at round-off (~1e-12).
            arguments
                obj
                eos   (1,1) thermo.ConfinedEOS
                Plist (:,1) double {mustBePositive}
                zlist (:,:) double {mustBeNonnegative}
            end
            if size(zlist, 1) ~= numel(obj.ComponentNames), zlist = zlist.'; end
            errLnPhi = 0;  errZ = 0;
            for P = Plist.'
                for j = 1:size(zlist, 2)
                    zj = zlist(:, j) / sum(zlist(:, j));
                    for flag = [-1 0 1]
                        [l1, Z1] = eos.calculateState(P, obj.T, zj, flag, obj.r_cap);
                        [l2, Z2] = obj.lnphi(P, zj, flag, eos.ShiftFugacity);
                        errLnPhi = max(errLnPhi, max(abs(l1 - l2)));
                        errZ = max(errZ, abs(Z1 - Z2));
                    end
                end
            end
        end

        function [lnphi, Z, V] = lnphi(obj, P, z, phaseFlag, shiftFugacity)
            % Independent standard-PR evaluation with the equivalent parameters
            % (no ConfinedEOS code path), used by verify().
            arguments
                obj
                P (1,1) double {mustBePositive}
                z (:,1) double
                phaseFlag (1,1) double {mustBeMember(phaseFlag, [-1 0 1])}
                shiftFugacity (1,1) logical = true
            end
            R_ = thermo.EquivalentPR.R;
            T_ = obj.T;
            w  = obj.omega(:);  Tc_ = obj.Tc(:);  Pc_ = obj.Pc(:);
            m = 0.37464 + 1.54226*w - 0.26992*w.^2;
            hi = w > 0.49;
            m(hi) = 0.379642 + 1.48503*w(hi) - 0.164423*w(hi).^2 + 0.016666*w(hi).^3;
            alpha = (1 + m .* (1 - sqrt(T_ ./ Tc_))).^2;
            a = obj.OmegaA(:) .* alpha .* (R_ * Tc_).^2 ./ Pc_;
            b = obj.OmegaB(:) .* R_ .* Tc_ ./ Pc_;
            Aij = sqrt(a * a.') .* (1 - obj.BIP);
            S = Aij * z;  am = z.' * S;  bm = z.' * b;
            RT = R_ * T_;
            A = am * P / RT^2;  B = bm * P / RT;
            zr = roots([1, B - 1, A - 3*B^2 - 2*B, B^3 + B^2 - A*B]);
            zr = sort(real(zr(abs(imag(zr)) <= 1e-10 * max(1, abs(real(zr))))));
            zr = zr(zr > B);
            d1 = 1 + sqrt(2);  d2 = 1 - sqrt(2);
            switch phaseFlag
                case  1, Z = zr(1);
                case -1, Z = zr(end);
                otherwise
                    g = zr - 1 - log(zr - B) - A/(2*sqrt(2)*B) * log((zr + d1*B) ./ (zr + d2*B));
                    [~, iMin] = min(g);
                    Z = zr(iMin);
            end
            lnphi = (b / bm) * (Z - 1) - log(Z - B) ...
                  - A / (2*sqrt(2)*B) * (2 * S / am - b / bm) * log((Z + d1*B) / (Z + d2*B));
            cs = obj.VolumeShift(:) .* b;
            V = Z * RT / P - z.' * cs;
            if shiftFugacity
                lnphi = lnphi - cs * P / RT;
            end
        end

        function writeTable(obj, filePath)
            % Plain per-component table plus the BIP' matrix (CSV), simulator-neutral.
            tbl = table(obj.ComponentNames.', obj.Tc.', obj.Pc.' / 1e5, obj.omega.', obj.MW.', ...
                obj.OmegaA.', obj.OmegaB.', obj.VolumeShift.', obj.CoverA.', ...
                'VariableNames', {'comp', 'Tc_K', 'Pc_bar', 'omega', 'MW', 'OmegaA', 'OmegaB', ...
                                  'SSHIFT', 'c_over_a'});
            writetable(tbl, filePath);
            [p, n] = fileparts(filePath);
            B = array2table(obj.BIP, 'VariableNames', cellstr(obj.ComponentNames), ...
                            'RowNames', cellstr(obj.ComponentNames));
            writetable(B, fullfile(p, n + "_BIP.csv"), 'WriteRowNames', true);
        end
    end

    methods (Static, Access = private)
        function checkCompatible(eqs)
            for k = 2:numel(eqs)
                if ~isequal(eqs(k).ComponentNames, eqs(1).ComponentNames)
                    error('EquivalentPR:ComponentMismatch', 'All EOS regions must share one component list.');
                end
                if abs(eqs(k).T - eqs(1).T) > 1e-9
                    error('EquivalentPR:TemperatureMismatch', 'All EOS regions must be built at the same T.');
                end
            end
        end

        function fid = openFile(filePath)
            d = fileparts(filePath);
            if strlength(d) > 0 && ~isfolder(d), mkdir(d); end
            fid = fopen(filePath, 'w');
            if fid < 0
                error('EquivalentPR:FileOpen', 'Cannot open %s for writing.', filePath);
            end
        end

        function e300Vec(fid, kw, eqs, getter)
            fprintf(fid, '\n%s\n', kw);
            for k = 1:numel(eqs)
                fprintf(fid, ' %.10g', getter(eqs(k)));
                fprintf(fid, ' /\n');
            end
        end

        function gemVec(fid, kw, v)
            fprintf(fid, '%s', kw);
            fprintf(fid, ' %.10g', v);
            fprintf(fid, '\n');
        end
    end
end
