classdef RockOffsetTIP
    % ROCKOFFSETTIP Frozen FWI-Only model: (k, p, lambda) and the rock-offset,
    % size-independent TIP correlation
    %
    %   ln(-TIP_ij) = alpha_rock + B * x_ij,      TIP_ij = -exp(...) < 0
    %
    % alpha_EF2, alpha_B1 and B are obtained by the same least-squares fit as
    % TIP_correlation.m (fitForm) from the frozen per-group optimised TIPs, with the
    % pair descriptors evaluated at the calibration temperature Tcal. At another
    % temperature the descriptor is re-evaluated at that T (only the Tstar family
    % depends on T); alpha and B stay frozen.
    %
    % Single source of the frozen values for PVT_verification, MLDO_check and
    % Temperature_test.
    %
    % Usage:
    %   tip = thermo.RockOffsetTIP.frozen("jL", workbookFile);
    %
    % Presets
    %   jL  : joint-Linearized, LL-inv_wbar (primary; workbook BIP C1-nC5 = 0.012022,
    %         adjusted C1-nC5 compositions)
    %   jL2 : jL with the C1-nC5 group TIP re-optimised for the recalibrated bulk
    %         model (nominal C1-nC5 compositions, BIP C1-nC5 = 0.0005, nC5-nC10 =
    %         0.000746; "MixtureData - Recalibrated.xlsx"). Same objective as the
    %         jL fit (RMSE of absolute confined Pdew, FWI-Only, k/p/lambda frozen);
    %         re-derived in MATLAB by Recalibrate_C1nC5.m.
    %   w1  : pareto-exact-w1, LL-Tstar (sensitivity; fails the temperature test)
    %   K   = tip.matrix(eos.Fluid, "EF2", T);           % TIP matrix (kijc)
    %   tip.applyTo(eos)                                  % set k, p, lambda on an EOS

    properties (Constant)
        kB   = 1.380649e-23        % [J/K]
        Tcal = 293.15              % [K] calibration temperature of the per-group TIPs
        Groups = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"]
    end

    properties (SetAccess = private)
        Preset      (1,1) string = ""
        Combination (1,1) string = ""
        k           (1,1) double = NaN
        p           (1,1) double = NaN
        lambda      (1,1) double = NaN
        Descriptor  (1,1) string = ""
        Form        (1,1) string = ""
        GroupTIPs   (1,:) double = zeros(1,0)   % per-group optimised TIPs (FWI-Only), order of Groups
        alphaEF2    (1,1) double = NaN
        alphaB1     (1,1) double = NaN
        B           (1,1) double = NaN
        rmseLnEF2   (1,1) double = NaN
        xCal        (:,1) double = zeros(0,1)
    end

    methods (Static)
        function obj = frozen(preset, workbookFile)
            % preset "jL" (joint-Linearized, LL-inv_wbar; primary) or
            % "w1" (pareto-exact-w1, LL-Tstar; sensitivity).
            arguments
                preset       (1,1) string {mustBeMember(preset, ["jL", "jL2", "w1"])}
                workbookFile (1,1) string
            end
            obj = thermo.RockOffsetTIP();
            obj.Preset = preset;
            switch preset
                case "jL"   % values as in TIP_correlation.m / Baclcalc_rp_B1_late.m
                    obj.Combination = "joint-linearized";
                    obj.k = 112.42634826677;  obj.p = 1.3558179290678;  obj.lambda = 2.0916235090532;
                    obj.Descriptor = "inv_wbar";  obj.Form = "LL-inv_wbar";
                    obj.GroupTIPs = [-4.9229696302796, -1.96920366521686, -1.44628069292089, -0.754862049795788];
                case "jL2"
                    obj.Combination = "joint-linearized (C1-nC5 recalibrated)";
                    obj.k = 112.42634826677;  obj.p = 1.3558179290678;  obj.lambda = 2.0916235090532;
                    obj.Descriptor = "inv_wbar";  obj.Form = "LL-inv_wbar";
                    obj.GroupTIPs = [-4.891452765367, -1.96920366521686, -1.44628069292089, -0.754862049795788];
                case "w1"
                    obj.Combination = "pareto-exact-w1";
                    obj.k = 166.410385934507;  obj.p = 1.5488421885245;  obj.lambda = 4.47908205631447;
                    obj.Descriptor = "Tstar";  obj.Form = "LL-Tstar";
                    obj.GroupTIPs = [-7.012375, -2.883946, -2.124570, -1.168401];
            end

            obj = thermo.RockOffsetTIP.fitCorrelation(obj, workbookFile);
        end

        function obj = custom(name, groupTIPs, descriptor, workbookFile, kpl)
            % Rock-offset correlation fitted to user-supplied per-group TIPs (order of
            % Groups), e.g. TIPs refitted at the edges of a BIP window. kpl = [k, p, lambda].
            arguments
                name         (1,1) string
                groupTIPs    (1,4) double {mustBeNegative}
                descriptor   (1,1) string
                workbookFile (1,1) string
                kpl          (1,3) double
            end
            obj = thermo.RockOffsetTIP();
            obj.Preset = "custom";  obj.Combination = name;
            obj.k = kpl(1);  obj.p = kpl(2);  obj.lambda = kpl(3);
            obj.Descriptor = descriptor;  obj.Form = "LL-" + descriptor;
            obj.GroupTIPs = groupTIPs;
            obj = thermo.RockOffsetTIP.fitCorrelation(obj, workbookFile);
        end

        function X = descriptor(fl, name, T)
            % Pair descriptors, identical to TIP_correlation.m pairDescriptor.
            w  = fl.omega(:);
            e  = fl.LJ_Energy(:);                       % [J]
            MW = fl.MW(:);
            Ts = (thermo.RockOffsetTIP.kB * T) ./ sqrt(e * e.');
            switch name
                case "wbar",      X = 0.5 * (w + w.');
                case "inv_wbar",  X = 2 ./ (w + w.');
                case "Tstar",     X = Ts;
                case "inv_Tstar", X = 1 ./ Ts;
                case "ln_Tstar",  X = log(Ts);
                case "inv_Mgeo",  X = 1 ./ sqrt(MW * MW.');
                otherwise, error('RockOffsetTIP:Descriptor', 'Unknown descriptor "%s".', name);
            end
        end
    end

    methods
        function K = matrix(obj, fl, rock, T)
            % TIP matrix for fluid fl on rock ("EF2" | "B1") at temperature T [K].
            switch string(rock)
                case "EF2", alpha = obj.alphaEF2;
                case "B1",  alpha = obj.alphaB1;
                otherwise,  error('RockOffsetTIP:Rock', 'No rock offset for "%s".', rock);
            end
            nc = fl.NC;
            K = -exp(alpha + obj.B * thermo.RockOffsetTIP.descriptor(fl, obj.Descriptor, T));
            K(1:nc+1:end) = 0;
            K = (K + K.') / 2;
        end

        function applyTo(obj, eos)
            % Set the frozen k, p, lambda on a ConfinedEOS (TIP matrix set separately).
            eos.k = obj.k;  eos.pT_wall = obj.p;  eos.lambda = obj.lambda;
            eos.ConfinementScale = 1;
        end

        function s = summary(obj)
            s = sprintf(['%s (%s): k = %.6f, p = %.6f, lambda = %.6f | TIP %s: alpha_EF2 = %+.6f, ' ...
                         'alpha_B1 = %+.6f, B = %+.6f (RMSE ln, EF2 = %.4f)'], obj.Preset, obj.Combination, ...
                         obj.k, obj.p, obj.lambda, obj.Form, obj.alphaEF2, obj.alphaB1, obj.B, obj.rmseLnEF2);
        end
    end

    methods (Static, Hidden)
        function obj = fitCorrelation(obj, workbookFile)
            % Least squares ln(-TIP_g) = alpha_rock(g) + B x_g over the four groups,
            % descriptors at Tcal (as TIP_correlation.m fitForm). Group descriptors
            % are cached per workbook and descriptor.
            persistent XC
            if isempty(XC), XC = containers.Map(); end
            G = numel(obj.Groups);
            ck = char(workbookFile + "|" + obj.Descriptor);
            if isKey(XC, ck)
                x = XC(ck);
            else
                x = zeros(G, 1);
                for g = 1:G
                    pr = split(obj.Groups(g), "|");
                    fl = entities.FluidProperties.loadFromWorkbook(workbookFile, pr(1));
                    X  = thermo.RockOffsetTIP.descriptor(fl, obj.Descriptor, obj.Tcal);
                    x(g) = X(1, 2);
                end
                XC(ck) = x;
            end
            isB1 = endsWith(obj.Groups(:), "|B1");
            y = log(-obj.GroupTIPs(:));
            A = [double(~isB1), double(isB1), x];
            if rank(A) < 3
                error('RockOffsetTIP:Rank', 'TIP design matrix is rank-deficient.');
            end
            par = A \ y;
            res = A * par - y;
            obj.alphaEF2 = par(1);  obj.alphaB1 = par(2);  obj.B = par(3);
            obj.rmseLnEF2 = sqrt(mean(res(~isB1).^2));
            obj.xCal = x;
        end
    end
end
