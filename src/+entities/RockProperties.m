classdef RockProperties
    % ROCKPROPERTIES Domain entity holding the per-mineral description of a
    % rock sample: wall energy parameter, grain density and mass fraction of
    % every XRD mineral and of the organic matter (TOC).
    %
    % Scope:
    %   - Reads, validates and stores raw per-mineral data. No grouping.
    %   - Exposes only rock-intrinsic basis transformations (mass closure,
    %     mass -> volume fraction). These depend on the rock alone.
    %   - Does NOT combine energies into an effective wall energy. That is a
    %     model (mixing-rule) decision and belongs to the EOS / wall model.
    %
    % Mass-basis convention of the Rock_Props sheet:
    %   XRD mineral fractions are reported on an inorganic (TOC-free) basis
    %   and sum to 1; TOC is reported on a bulk-rock basis. The bulk
    %   composition is therefore
    %       w_org   = f_k * TOC
    %       w_i     = w_i,XRD / sum(w_XRD) * (1 - w_org)
    %   where f_k is the kerogen/TOC mass conversion factor (default 1.0,
    %   i.e. TOC treated as the organic mass). Data already on a bulk basis
    %   (all entries sum to 1) map to the same expression.

    properties (SetAccess = private)
        RockName        (1,1) string = ""
        Theta           (1,1) double {mustBeReal, mustBeNonnegative} = 0   % contact angle [deg]

        MineralNames    (1,:) string = string.empty(1,0)  % includes "TOC"
        Epsilon_K       (1,:) double = zeros(1,0)          % energy parameter [K]
        GrainDensity    (1,:) double = zeros(1,0)          % [g/cm^3]
        ReportedMassFraction (1,:) double = zeros(1,0)     % as read from sheet

        RawTOC          (1,1) double {mustBeReal, mustBeNonnegative} = 0   % bulk TOC mass fraction (TIP correlation input)
        KerogenFactor   (1,1) double {mustBeReal, mustBePositive} = 1.0   % kerogen mass / TOC mass
        ReportedBasis   (1,1) string = ""                  % "inorganic" | "bulk"
    end

    properties (Constant)
        R = 8.3144621                   % [J/(mol K)]
        OrganicLabel = "TOC"
    end

    properties (Dependent)
        NumMinerals                     % scalar
        IsOrganic   % mask of the TOC column
        Energy_Jmol    % Epsilon_K * R
        MassFraction    % bulk basis, sums to 1
        VolumeFraction    % bulk basis, sums to 1
    end

    methods
        function obj = RockProperties(name, theta, mineralNames, epsilonK, rho, wReported, opts)
            arguments
                name = ""
                theta (1,1) double = 0
                mineralNames (1,:) string = string.empty(1,0)
                epsilonK (1,:) double = zeros(1,0)
                rho (1,:) double = zeros(1,0)
                wReported (1,:) double = zeros(1,0)
                opts.KerogenFactor (1,1) double {mustBePositive} = 1.0
                opts.ClosureTol (1,1) double {mustBePositive} = 5e-3
            end
            if nargin == 0, return; end

            n = numel(mineralNames);
            if any([numel(epsilonK), numel(rho), numel(wReported)] ~= n)
                error('RockProperties:SizeMismatch', ...
                    'Energy, density and mass-fraction vectors must match the %d mineral columns.', n);
            end

            % Absent minerals: blank weight -> 0; blank energy/density only
            % tolerated if the mineral is absent from this rock.
            wReported(isnan(wReported)) = 0;
            present = wReported > 0;
            if any(present & (isnan(epsilonK) | isnan(rho)))
                bad = strjoin(mineralNames(present & (isnan(epsilonK) | isnan(rho))), ', ');
                error('RockProperties:MissingParameter', ...
                    'Energy or grain density missing for present mineral(s): %s', bad);
            end
            if any(wReported < 0) || any(epsilonK(present) < 0) || any(rho(present) <= 0)
                error('RockProperties:InvalidValue', 'Negative weight/energy or non-positive density.');
            end

            isOrg = mineralNames == entities.RockProperties.OrganicLabel;
            if nnz(isOrg) > 1
                error('RockProperties:DuplicateTOC', 'More than one "%s" column.', entities.RockProperties.OrganicLabel);
            end
            toc = sum(wReported(isOrg));
            sInorg = sum(wReported(~isOrg));

            % Identify reporting basis (validation + traceability only)
            if abs(sInorg - 1) <= opts.ClosureTol
                basis = "inorganic";
            elseif abs(sInorg + toc - 1) <= opts.ClosureTol
                basis = "bulk";
            else
                error('RockProperties:ClosureFailure', ...
                    ['%s: inorganic sum = %.4f, total = %.4f. Neither the TOC-free ' ...
                     'nor the bulk basis closes within %.1e.'], string(name), sInorg, sInorg + toc, opts.ClosureTol);
            end
            if opts.KerogenFactor * toc >= 1
                error('RockProperties:OrganicOverflow', 'Kerogen mass fraction >= 1.');
            end

            obj.RockName      = string(name);
            obj.Theta         = theta;
            obj.MineralNames  = mineralNames;
            obj.Epsilon_K     = epsilonK;
            obj.GrainDensity  = rho;
            obj.ReportedMassFraction = wReported;
            obj.RawTOC        = toc;
            obj.KerogenFactor = opts.KerogenFactor;
            obj.ReportedBasis = basis;
        end

        % ---------------- Dependent getters ----------------
        function n = get.NumMinerals(obj)
            n = numel(obj.MineralNames);
        end

        function m = get.IsOrganic(obj)
            m = obj.MineralNames == entities.RockProperties.OrganicLabel;
        end

        function E = get.Energy_Jmol(obj)
            E = obj.Epsilon_K * entities.RockProperties.R;
        end

        function w = get.MassFraction(obj)
            w = zeros(1, obj.NumMinerals);
            if obj.NumMinerals == 0, return; end
            org = obj.IsOrganic;
            wOrg = obj.KerogenFactor * obj.RawTOC;
            wIn  = obj.ReportedMassFraction(~org);
            w(~org) = wIn / sum(wIn) * (1 - wOrg);
            w(org)  = wOrg;
        end

        function v = get.VolumeFraction(obj)
            w = obj.MassFraction;
            v = zeros(size(w));
            p = w > 0;                       % avoid NaN density of absent minerals
            v(p) = w(p) ./ obj.GrainDensity(p);
            v = v / sum(v);
        end

        % ---------------- Accessors ----------------
        function s = mineral(obj, name)
            % Returns struct of all per-mineral data for one mineral.
            idx = find(obj.MineralNames == string(name), 1);
            if isempty(idx)
                error('RockProperties:UnknownMineral', '"%s" not in %s.', string(name), obj.RockName);
            end
            s = struct('Name', obj.MineralNames(idx), ...
                       'Epsilon_K', obj.Epsilon_K(idx), ...
                       'Energy_Jmol', obj.Energy_Jmol(idx), ...
                       'GrainDensity', obj.GrainDensity(idx), ...
                       'MassFraction', obj.MassFraction(idx), ...
                       'VolumeFraction', obj.VolumeFraction(idx));
        end

        function T = toTable(obj)
            % Diagnostic view of the rock composition.
            T = table(obj.MineralNames.', obj.Epsilon_K.', obj.Energy_Jmol.', ...
                      obj.GrainDensity.', obj.ReportedMassFraction.', ...
                      obj.MassFraction.', obj.VolumeFraction.', ...
                'VariableNames', {'Mineral','Epsilon_K','E_Jmol','Rho_gcc', ...
                                  'w_reported','w_bulk','v_bulk'});
        end
    end

    methods (Static)
        function obj = loadFromWorkbook(filePath, targetRockName, opts)
            % Reads the Rock_Props sheet. Mineral columns are discovered from
            % the header (every column except Mineral / Rock / theta), so
            % adding a mineral requires no code change.
            arguments
                filePath (1,1) string
                targetRockName (1,1) string
                opts.Sheet (1,1) string = "Rock_Props"
                opts.KerogenFactor (1,1) double {mustBePositive} = 1.0
                opts.ClosureTol (1,1) double {mustBePositive} = 5e-3
            end

            io = detectImportOptions(filePath, 'Sheet', opts.Sheet);
            io.VariableNamingRule = 'preserve';
            T = readtable(filePath, io);

            vn = string(T.Properties.VariableNames);
            required = ["Mineral", "Rock", "theta"];
            if ~all(ismember(required, vn))
                error('RockProperties:InvalidSheetFormat', ...
                    'Sheet "%s" must contain columns: %s.', opts.Sheet, strjoin(required, ', '));
            end
            mineralNames = vn(~ismember(vn, required));

            rowType = strtrim(string(T.("Mineral")));
            rockCol = strtrim(string(T.("Rock")));

            iE   = entities.RockProperties.uniqueRow(rowType == "mineral_E",   "mineral_E");
            iRho = entities.RockProperties.uniqueRow(rowType == "mineral_rho", "mineral_rho");
            iW   = entities.RockProperties.uniqueRow(rowType == "mineral_w" & rockCol == targetRockName, ...
                                             "mineral_w for rock " + targetRockName);

            epsK = entities.RockProperties.rowValues(T, iE,   mineralNames);
            rho  = entities.RockProperties.rowValues(T, iRho, mineralNames);
            w    = entities.RockProperties.rowValues(T, iW,   mineralNames);
            theta = entities.RockProperties.toDouble(T.("theta")(iW));

            obj = entities.RockProperties(targetRockName, theta, mineralNames, epsK, rho, w, ...
                                 'KerogenFactor', opts.KerogenFactor, ...
                                 'ClosureTol', opts.ClosureTol);
        end
    end

    methods (Static, Access = private)
        function idx = uniqueRow(mask, label)
            idx = find(mask);
            if isempty(idx)
                error('RockProperties:RowNotFound', 'No "%s" row in sheet.', label);
            elseif numel(idx) > 1
                error('RockProperties:DuplicateRow', 'Multiple "%s" rows in sheet.', label);
            end
        end

        function vals = rowValues(T, idx, cols)
            vals = zeros(1, numel(cols));
            for k = 1:numel(cols)
                vals(k) = entities.RockProperties.toDouble(T.(cols(k))(idx));
            end
        end

        function x = toDouble(v)
            if iscell(v), v = v{1}; end
            if isnumeric(v) || islogical(v)
                x = double(v);
            else
                x = str2double(string(v));
            end
            if isempty(x) || ismissing(x), x = NaN; end
        end
    end
end
