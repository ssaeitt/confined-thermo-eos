classdef DewPointDataset
    % DEWPOINTDATASET Experimental bulk and confined dew points (Mix_Props).
    % Pressures are stored as in the sheet [psi]; getCase returns SI values.
    % The confined value is an interval [lower, upper]; its midpoint is the
    % target and its half-width the experimental uncertainty.

    properties (Constant)
        psi2Pa = 6894.757
    end

    properties (SetAccess = private)
        Mixture   (:,1) string = string.empty(0,1)
        MoleFrac  (:,1) cell   = cell(0,1)        % each a column vector
        Rock      (:,1) string = string.empty(0,1)
        T         (:,1) double = zeros(0,1)        % [K]
        BulkPdew_psi   (:,1) double = zeros(0,1)
        ConfLower_psi  (:,1) double = zeros(0,1)   % NaN if not measured
        ConfUpper_psi  (:,1) double = zeros(0,1)
        PoreRadius_nm  (:,1) double = zeros(0,1)   % optional column; NaN if absent
        Campaign       (:,1) string = string.empty(0,1)   % optional column; "" if absent
    end

    properties (Dependent)
        NumCases
        NumComponents       % per case
        HasConfined         % logical per case
        ConfMid_psi
        ConfHalfWidth_psi
        ShiftMid_psi        % confined midpoint - bulk
    end

    methods
        function obj = DewPointDataset(mix, z, rock, T, bulk, lo, up, rp_nm, campaign)
            if nargin == 0, return; end
            if nargin < 8 || isempty(rp_nm),    rp_nm = nan(numel(mix), 1); end
            if nargin < 9 || isempty(campaign), campaign = strings(numel(mix), 1); end
            obj.PoreRadius_nm = rp_nm(:);
            obj.Campaign = string(campaign(:));
            obj.Mixture = string(mix(:));
            obj.MoleFrac = z(:);
            obj.Rock = string(rock(:));
            obj.T = T(:);
            obj.BulkPdew_psi = bulk(:);
            obj.ConfLower_psi = lo(:);
            obj.ConfUpper_psi = up(:);
            for i = 1:numel(obj.Mixture)
                nComp = numel(split(obj.Mixture(i), "-"));
                zi = obj.MoleFrac{i};
                if numel(zi) ~= nComp
                    error('DewPointDataset:CompositionSize', ...
                        'Row %d (%s): %d mole fractions for %d components.', i, obj.Mixture(i), numel(zi), nComp);
                end
                if any(~isfinite(zi)) || abs(sum(zi) - 1) > 1e-6
                    error('DewPointDataset:CompositionSum', ...
                        'Row %d (%s): mole fractions sum to %.6f.', i, obj.Mixture(i), sum(zi));
                end
            end
        end

        function n = get.NumCases(obj),       n = numel(obj.Mixture); end
        function n = get.NumComponents(obj),  n = cellfun(@numel, obj.MoleFrac); end
        function m = get.HasConfined(obj),    m = isfinite(obj.ConfLower_psi) & isfinite(obj.ConfUpper_psi); end
        function v = get.ConfMid_psi(obj),    v = 0.5 * (obj.ConfLower_psi + obj.ConfUpper_psi); end
        function v = get.ConfHalfWidth_psi(obj), v = 0.5 * (obj.ConfUpper_psi - obj.ConfLower_psi); end
        function v = get.ShiftMid_psi(obj),   v = obj.ConfMid_psi - obj.BulkPdew_psi; end

        function c = getCase(obj, i)
            % Single case in SI units.
            p = entities.DewPointDataset.psi2Pa;
            c = struct( ...
                'index',     i, ...
                'Mixture',   obj.Mixture(i), ...
                'Components', split(obj.Mixture(i), "-").', ...
                'z',         obj.MoleFrac{i}, ...
                'Rock',      obj.Rock(i), ...
                'T',         obj.T(i), ...
                'PdewBulk',  obj.BulkPdew_psi(i) * p, ...
                'PdewConfMid', obj.ConfMid_psi(i) * p, ...
                'PdewConfHalfWidth', obj.ConfHalfWidth_psi(i) * p, ...
                'HasConfined', obj.HasConfined(i), ...
                'PoreRadius_m', obj.PoreRadius_nm(i) * 1e-9, ...   % NaN when the sheet has no column
                'Campaign',   obj.Campaign(i));
        end

        function idx = findCases(obj, mixture, rock, z)
            % Row indices matching mixture/rock and, optionally, composition.
            idx = find(obj.Mixture == string(mixture) & obj.Rock == string(rock));
            if nargin > 3 && ~isempty(z)
                keep = arrayfun(@(k) numel(obj.MoleFrac{k}) == numel(z) && ...
                    max(abs(obj.MoleFrac{k} - z(:))) < 1e-9, idx);
                idx = idx(keep);
            end
        end

        function Tb = toTable(obj)
            zs = cellfun(@(v) "[" + strjoin(string(v.'), ",") + "]", obj.MoleFrac);
            Tb = table(obj.Mixture, zs, obj.Rock, obj.T, obj.BulkPdew_psi, ...
                obj.ConfLower_psi, obj.ConfUpper_psi, obj.PoreRadius_nm, obj.Campaign, ...
                'VariableNames', {'Mixture','z','Rock','T_K','Bulk_psi','ConfLo_psi','ConfUp_psi', ...
                'rp_nm','Campaign'});
        end
    end

    methods (Static)
        function obj = loadFromWorkbook(filePath, opts)
            arguments
                filePath (1,1) string
                opts.Sheet (1,1) string = "Mix_Props"
                opts.DefaultT (1,1) double = 293.15    % used when the sheet has no T column
            end
            io = detectImportOptions(filePath, 'Sheet', opts.Sheet);
            io.VariableNamingRule = 'preserve';
            Tb = readtable(filePath, io);
            vn = string(Tb.Properties.VariableNames);
            need = ["Mixture", "mole_frac", "Bulk_Pdew", "Rock", "Pdew_Exp_lower", "Pdew_Exp_upper"];
            if ~all(ismember(need, vn))
                error('DewPointDataset:InvalidSheetFormat', 'Sheet "%s" must contain: %s.', ...
                    opts.Sheet, strjoin(need, ', '));
            end

            mix  = strtrim(string(Tb.("Mixture")));
            keep = ~ismissing(mix) & mix ~= "";
            Tb = Tb(keep, :);
            mix = mix(keep);

            zStr = string(Tb.("mole_frac"));
            z = cell(numel(mix), 1);
            for i = 1:numel(mix)
                z{i} = str2double(split(erase(zStr(i), ["[", "]", " "]), ","));
            end

            if ismember("T", vn)
                T = entities.DewPointDataset.num(Tb.("T"));
            else
                T = repmat(opts.DefaultT, numel(mix), 1);
            end

            % Optional columns: pore radius [nm] and campaign/batch label
            rp = nan(numel(mix), 1);
            rpCol = vn(ismember(lower(vn), ["r_p_nm", "rp_nm", "pore_radius_nm", "rp"]));
            if ~isempty(rpCol), rp = entities.DewPointDataset.num(Tb.(rpCol(1))); end
            camp = strings(numel(mix), 1);
            cCol = vn(ismember(lower(vn), ["campaign", "batch", "note"]));
            if ~isempty(cCol), camp = strtrim(string(Tb.(cCol(1)))); end

            obj = entities.DewPointDataset(mix, z, strtrim(string(Tb.("Rock"))), T, ...
                entities.DewPointDataset.num(Tb.("Bulk_Pdew")), ...
                entities.DewPointDataset.num(Tb.("Pdew_Exp_lower")), ...
                entities.DewPointDataset.num(Tb.("Pdew_Exp_upper")), rp, camp);
        end
    end

    methods (Static, Access = private)
        function x = num(v)
            if isnumeric(v), x = double(v); else, x = str2double(string(v)); end
            x = x(:);
        end
    end
end
