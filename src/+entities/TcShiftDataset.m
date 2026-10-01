classdef TcShiftDataset
    % TCSHIFTDATASET Literature confined critical-point shifts of pure
    % components (calibration data for the FWI parameters k, lambda, pT_wall).
    %
    %   dTc = (Tc,bulk - Tc,pore) / Tc,bulk,   dPc likewise (not fitted: the
    %   model has an unconfined b, so it predicts dPc = dTc).
    %
    % Geometry convention. ConfinedEOS uses dH = sigma_i / r_cap with r_cap the
    % pore RADIUS. The tabulated d_H values are not homogeneous:
    %   Kuang et al (2025)          sigma / R_p  (cylinder radius)  -> factor 1
    %   Pitakbunkate et al (2016)   sigma / H    (slit width)       -> factor 2 for a radius basis
    %   Singh et al (2009)          sigma / H    (slit width)       -> factor 2
    %   Sobecki et al (2019)        sigma / H    (slit width)       -> factor 2
    % GeometryFactor (per reference, default 1 = tabulated values used as-is)
    % converts tabulated d_H to the model's sigma/r basis: dH_model = factor * d_H.

    properties (SetAccess = private)
        Component  (:,1) string = string.empty(0,1)
        dH         (:,1) double = zeros(0,1)     % as tabulated
        dTc        (:,1) double = zeros(0,1)
        dPc        (:,1) double = zeros(0,1)     % NaN where not reported
        Reference  (:,1) string = string.empty(0,1)
        GeometryFactor (:,1) double = zeros(0,1) % per row
        SourceFile  (1,1) string = ""
        SourceSheet (1,1) string = ""
    end

    properties (Dependent)
        NumPoints
        dH_model                                 % GeometryFactor .* dH
    end

    methods
        function obj = TcShiftDataset(component, dH, dTc, dPc, reference, factor, src, sheet)
            if nargin == 0, return; end
            n = numel(component);
            if any([numel(dH), numel(dTc), numel(dPc), numel(reference), numel(factor)] ~= n)
                error('TcShiftDataset:SizeMismatch', 'All columns must have the same length.');
            end
            if any(~isfinite(dH) | dH <= 0) || any(~isfinite(dTc))
                error('TcShiftDataset:InvalidValue', 'd_H must be positive and dTc finite.');
            end
            obj.Component = strtrim(string(component(:)));
            obj.dH  = dH(:);
            obj.dTc = dTc(:);
            obj.dPc = dPc(:);
            obj.Reference = strtrim(string(reference(:)));
            obj.GeometryFactor = factor(:);
            obj.SourceFile = string(src);
            obj.SourceSheet = string(sheet);
        end

        function n = get.NumPoints(obj)
            n = numel(obj.dTc);
        end

        function v = get.dH_model(obj)
            v = obj.GeometryFactor .* obj.dH;
        end

        function sub = select(obj, mask)
            sub = entities.TcShiftDataset(obj.Component(mask), obj.dH(mask), obj.dTc(mask), ...
                obj.dPc(mask), obj.Reference(mask), obj.GeometryFactor(mask), ...
                obj.SourceFile, obj.SourceSheet);
        end

        function [compIdx, r_cap] = modelInputs(obj, fluid)
            % Component index in FluidProperties and the pore radius [m] that
            % reproduces dH_model with the model's own sigma: r = sigma_i / dH_model.
            compIdx = zeros(obj.NumPoints, 1);
            for n = 1:obj.NumPoints
                compIdx(n) = fluid.componentIndex(obj.Component(n));
            end
            sigma = fluid.LJ_Size(:);
            r_cap = sigma(compIdx) ./ obj.dH_model;
        end

        function T = toTable(obj)
            T = table(obj.Component, obj.dH, obj.dH_model, obj.dTc, obj.dPc, obj.Reference, ...
                'VariableNames', {'Component','dH','dH_model','dTc','dPc','Reference'});
        end
    end

    methods (Static)
        function obj = loadFromWorkbook(filePath, opts)
            arguments
                filePath (1,1) string
                opts.Sheet (1,1) string = "All Data"
                opts.GeometryReferences (1,:) string = string.empty(1,0)
                opts.GeometryFactors    (1,:) double = zeros(1,0)
                opts.Components         (1,:) string = string.empty(1,0)  % empty = all
            end
            if numel(opts.GeometryReferences) ~= numel(opts.GeometryFactors)
                error('TcShiftDataset:GeometrySpec', ...
                    'GeometryReferences and GeometryFactors must have equal length.');
            end

            io = detectImportOptions(filePath, 'Sheet', opts.Sheet);
            io.VariableNamingRule = 'preserve';
            T = readtable(filePath, io);

            need = ["Component", "d_H", "delta_Tc", "delta_Pc", "Reference"];
            vn = string(T.Properties.VariableNames);
            if ~all(ismember(need, vn))
                error('TcShiftDataset:InvalidSheetFormat', 'Sheet "%s" must contain: %s.', ...
                    opts.Sheet, strjoin(need, ', '));
            end

            comp = strtrim(string(T.("Component")));
            dH   = entities.TcShiftDataset.num(T.("d_H"));
            dTc  = entities.TcShiftDataset.num(T.("delta_Tc"));
            dPc  = entities.TcShiftDataset.num(T.("delta_Pc"));
            ref  = strtrim(string(T.("Reference")));

            keep = ~ismissing(comp) & comp ~= "" & isfinite(dH) & isfinite(dTc);
            if ~isempty(opts.Components)
                keep = keep & ismember(comp, opts.Components);
            end

            factor = ones(size(dH));
            for k = 1:numel(opts.GeometryReferences)
                hit = ref == opts.GeometryReferences(k);
                if ~any(hit & keep)
                    warning('TcShiftDataset:UnusedGeometryFactor', ...
                        'No rows for reference "%s".', opts.GeometryReferences(k));
                end
                factor(hit) = opts.GeometryFactors(k);
            end

            obj = entities.TcShiftDataset(comp(keep), dH(keep), dTc(keep), dPc(keep), ...
                ref(keep), factor(keep), filePath, opts.Sheet);
        end
    end

    methods (Static, Access = private)
        function x = num(v)
            if isnumeric(v)
                x = double(v);
            else
                x = str2double(string(v));
            end
            x = x(:);
        end
    end
end
