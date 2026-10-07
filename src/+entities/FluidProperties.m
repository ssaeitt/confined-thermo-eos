classdef FluidProperties
    % FLUIDPROPERTIES Domain entity encapsulating pure-component constants,
    % volume-translation parameters, Lennard-Jones parameters and the
    % fluid-fluid binary interaction matrix (BIP) for an NC-component mixture.
    %
    % Sheet conventions (MixtureData.xlsx):
    %   Comp_Props : one row per component, header row with
    %                comp | Tc [K] | Pc [Pa] | omega | Mw [g/mol] | Vc | parachor |
    %                Zrackett | SE | LJ_Size [m] | LJ_Energy [J]
    %   BIP        : square table, row labels in column 1, column labels in
    %                header row. Rows and columns are mapped independently by
    %                name, so their ordering in the sheet is irrelevant.

    properties (SetAccess = private)
        ComponentNames (1,:) string = string.empty(1,0)
        NC (1,1) double {mustBeInteger, mustBeNonnegative} = 0

        % Core cubic EOS parameters
        Tc    (1,:) double = zeros(1,0)   % [K]
        Pc    (1,:) double = zeros(1,0)   % [Pa]
        omega (1,:) double = zeros(1,0)   % [-]
        MW    (1,:) double = zeros(1,0)   % [g/mol]
        Vc    (1,:) double = zeros(1,0)   % as stored in sheet (see unit note in loader)

        % Liquid volume and interfacial tension
        Parachor    (1,:) double = zeros(1,0)   % MacLeod-Sugden parachor
        ZRackett    (1,:) double = zeros(1,0)   % Rackett compressibility
        VolumeShift (1,:) double = zeros(1,0)   % Peneloux shift, dimensionless s_i = c_i/b_i

        % Lennard-Jones parameters (fluid-wall interaction)
        LJ_Size   (1,:) double = zeros(1,0)     % sigma [m]
        LJ_Energy (1,:) double = zeros(1,0)     % epsilon [J]

        % Fluid-fluid cross interactions
        BIPMatrix (:,:) double = zeros(0,0)     % symmetric, zero diagonal [NC x NC]
    end

    properties (Constant, Access = private)
        kB = 1.380649e-23                       % [J/K]
        SymTol = 1e-12
    end

    properties (Dependent)
        LJ_Energy_K                             % epsilon/kB [K], diagnostic
    end

    methods
        function obj = FluidProperties(names, tc, pc, omega, mw, vc, parachor, zra, vshift, lj_size, lj_energy, bip)
            if nargin == 0, return; end

            obj.ComponentNames = reshape(strtrim(string(names)), 1, []);
            obj.NC = numel(obj.ComponentNames);
            if numel(unique(obj.ComponentNames)) ~= obj.NC
                error('FluidProperties:DuplicateComponent', 'Duplicate component names.');
            end

            obj.Tc          = reshape(double(tc), 1, []);
            obj.Pc          = reshape(double(pc), 1, []);
            obj.omega       = reshape(double(omega), 1, []);
            obj.MW          = reshape(double(mw), 1, []);
            obj.Vc          = reshape(double(vc), 1, []);
            obj.Parachor    = reshape(double(parachor), 1, []);
            obj.ZRackett    = reshape(double(zra), 1, []);
            obj.VolumeShift = reshape(double(vshift), 1, []);
            obj.LJ_Size     = reshape(double(lj_size), 1, []);
            obj.LJ_Energy   = reshape(double(lj_energy), 1, []);
            obj.BIPMatrix   = double(bip);

            obj.validate();
        end

        function e = get.LJ_Energy_K(obj)
            e = obj.LJ_Energy / entities.FluidProperties.kB;
        end

        function idx = componentIndex(obj, name)
            idx = find(obj.ComponentNames == string(name), 1);
            if isempty(idx)
                error('FluidProperties:UnknownComponent', '"%s" not in mixture.', string(name));
            end
        end

        function T = toTable(obj)
            % Diagnostic view.
            T = table(obj.ComponentNames.', obj.Tc.', obj.Pc.', obj.omega.', obj.MW.', ...
                      obj.Vc.', obj.Parachor.', obj.ZRackett.', obj.VolumeShift.', ...
                      obj.LJ_Size.', obj.LJ_Energy_K.', ...
                'VariableNames', {'comp','Tc_K','Pc_Pa','omega','Mw','Vc','parachor', ...
                                  'Zrackett','SE','sigma_m','epsK_K'});
        end
    end

    methods (Access = private)
        function validate(obj)
            n = obj.NC;
            vecs  = {obj.Tc, obj.Pc, obj.omega, obj.MW, obj.Vc, obj.Parachor, ...
                     obj.ZRackett, obj.VolumeShift, obj.LJ_Size, obj.LJ_Energy};
            names = ["Tc","Pc","omega","Mw","Vc","parachor","Zrackett","SE","LJ_Size","LJ_Energy"];

            for k = 1:numel(vecs)
                v = vecs{k};
                if numel(v) ~= n
                    error('FluidProperties:VectorSizingFailure', ...
                        '"%s" has %d entries; expected NC = %d.', names(k), numel(v), n);
                end
                bad = ~isfinite(v);
                if any(bad)
                    error('FluidProperties:MissingData', '"%s" missing/non-finite for: %s.', ...
                        names(k), strjoin(obj.ComponentNames(bad), ', '));
                end
            end

            % Sign constraints (omega and SE may be negative)
            positive = {obj.Tc, obj.Pc, obj.MW, obj.Vc, obj.Parachor, obj.ZRackett, ...
                        obj.LJ_Size, obj.LJ_Energy};
            pNames   = ["Tc","Pc","Mw","Vc","parachor","Zrackett","LJ_Size","LJ_Energy"];
            for k = 1:numel(positive)
                if any(positive{k} <= 0)
                    error('FluidProperties:NonPositive', '"%s" must be strictly positive.', pNames(k));
                end
            end

            % BIP: size, finiteness, symmetry, zero diagonal
            K = obj.BIPMatrix;
            if ~isequal(size(K), [n n])
                error('FluidProperties:BIPDimensionMismatch', ...
                    'BIP matrix is %dx%d; expected %dx%d.', size(K,1), size(K,2), n, n);
            end
            if any(~isfinite(K), 'all')
                error('FluidProperties:BIPMissing', 'BIP matrix contains missing/non-finite entries.');
            end
            if max(abs(K - K.'), [], 'all') > entities.FluidProperties.SymTol
                error('FluidProperties:BIPAsymmetric', 'BIP matrix is not symmetric.');
            end
            if any(abs(diag(K)) > entities.FluidProperties.SymTol)
                error('FluidProperties:BIPDiagonal', 'BIP diagonal must be zero.');
            end
        end
    end

    methods (Static)
        function obj = loadFromWorkbook(filePath, selection)
            % selection: mixture name "C1-nC5-nC10" (order preserved, matches the
            % mole-fraction vector in Mix_Props), a string array of component
            % names, or omitted/empty to load every component in Comp_Props.
            arguments
                filePath (1,1) string
                selection string = string.empty
            end

            % ---------- 1. Pure-component properties ----------
            opts = detectImportOptions(filePath, 'Sheet', 'Comp_Props');
            opts.VariableNamingRule = 'preserve';
            opts = setvartype(opts, 'comp', 'string');
            P = readtable(filePath, opts);

            required = ["comp","Tc","Pc","omega","Mw","Vc","parachor","Zrackett", ...
                        "SE","LJ_Size","LJ_Energy"];
            missingCols = required(~ismember(required, string(P.Properties.VariableNames)));
            if ~isempty(missingCols)
                error('FluidProperties:InvalidSheetFormat', ...
                    'Comp_Props is missing column(s): %s.', strjoin(missingCols, ', '));
            end

            % Drop blank/formatted-but-empty rows
            compCol = strtrim(P.("comp"));
            P = P(~ismissing(compCol) & compCol ~= "", :);
            compCol = strtrim(P.("comp"));

            % Component selection (explicit order)
            if isempty(selection) || all(selection == "")
                target = reshape(compCol, 1, []);
            elseif isscalar(selection)
                target = strtrim(split(selection, "-")).';
            else
                target = reshape(strtrim(selection), 1, []);
            end

            rowIdx = zeros(1, numel(target));
            for i = 1:numel(target)
                hit = find(compCol == target(i));
                if isempty(hit)
                    error('FluidProperties:ComponentNotFound', ...
                        'Component "%s" not found in Comp_Props.', target(i));
                elseif numel(hit) > 1
                    error('FluidProperties:DuplicateComponent', ...
                        'Component "%s" appears %d times in Comp_Props.', target(i), numel(hit));
                end
                rowIdx(i) = hit;
            end
            P = P(rowIdx, :);

            % ---------- 2. BIP matrix (rows and columns mapped by name) ----------
            optsB = detectImportOptions(filePath, 'Sheet', 'BIP');
            optsB.VariableNamingRule = 'preserve';
            B = readtable(filePath, optsB);

            rowLabels = strtrim(string(B{:, 1}));
            colLabels = strtrim(string(B.Properties.VariableNames(2:end)));
            valid = ~ismissing(rowLabels) & rowLabels ~= "";
            B = B(valid, :);
            rowLabels = rowLabels(valid);
            Kraw = double(B{:, 2:end});

            [okR, rIdx] = ismember(target, rowLabels);
            [okC, cIdx] = ismember(target, colLabels);
            if ~all(okR) || ~all(okC)
                miss = unique([target(~okR), target(~okC)]);
                error('FluidProperties:BIPComponentNotFound', ...
                    'Component(s) absent from BIP row/column labels: %s.', strjoin(miss, ', '));
            end
            K = Kraw(rIdx, cIdx);

            % ---------- 3. Instantiate ----------
            obj = entities.FluidProperties(target, P.("Tc"), P.("Pc"), P.("omega"), ...
                P.("Mw"), P.("Vc"), P.("parachor"), P.("Zrackett"), P.("SE"), ...
                P.("LJ_Size"), P.("LJ_Energy"), K);
        end
    end
end
