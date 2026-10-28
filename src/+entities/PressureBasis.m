classdef PressureBasis
    % PRESSUREBASIS Gauge / absolute basis of the recorded pressures.
    %
    % The workbook keeps every pressure as recorded (Doryani thesis tables, late
    % campaign values). A per-row flag column states the basis of the record:
    %
    %   P_basis = "psig"  gauge, early campaign (LabView VI / Druck display unit)
    %             "psia"  absolute, late campaign (Quartzdyne with Q-Com, from ~2017)
    %
    % On loading with Mode "absolute" the gauge rows get + Atm_psi on every column in
    % Columns (bulk Pdew, confined range, P at MLDO). Shifts (confined - bulk) are
    % differences of two values on the same basis and are unchanged. Mode "recorded"
    % leaves the values as recorded (reproduces the results obtained before the
    % basis was known).
    %
    % Atm_psi is the standard atmosphere (101.325 kPa). The local barometric
    % pressure during the tests is not recorded; its spread (about +-0.3 psi over a
    % year at sea level) is an uncertainty of the conversion, small against the
    % other terms of the budget.
    %
    % Usage:
    %   T  = entities.PressureBasis.readSheet(file, "Bulk_PVT", "absolute");
    %   T  = entities.PressureBasis.apply(T, "absolute", "Bulk_PVT");
    %   off = entities.PressureBasis.offset(T, "absolute", "Mix_Props");   % [psi] per row
    %
    % tools/add_pressure_basis.m writes the P_basis column into a workbook.

    properties (Constant)
        Atm_psi     = 14.695949        % 101.325 kPa
        FlagColumn  = "P_basis"
        Columns     = ["Bulk_Pdew", "Pdew_Exp_lower", "Pdew_Exp_upper", "P_MLDO_psi"]
        Modes       = ["absolute", "recorded"]
    end

    methods (Static)
        function off = offset(Tb, mode, sheet)
            % Offset [psi] to add to every pressure column of each row.
            arguments
                Tb    table
                mode  (1,1) string {mustBeMember(mode, ["absolute", "recorded"])}
                sheet (1,1) string = "(table)"
            end
            off = zeros(height(Tb), 1);
            if mode == "recorded", return; end
            vn = string(Tb.Properties.VariableNames);
            if ~ismember(entities.PressureBasis.FlagColumn, vn)
                error('PressureBasis:NoFlag', ['Sheet "%s" has no %s column; absolute pressures need it. ' ...
                    'Run tools/add_pressure_basis.m on the workbook, or load with mode "recorded".'], ...
                    sheet, entities.PressureBasis.FlagColumn);
            end
            b = lower(strtrim(string(Tb.(entities.PressureBasis.FlagColumn))));
            b(ismissing(b)) = "";
            % Blank rows (e.g. padding below the data) are dropped by the loaders;
            % any other value must be psig or psia.
            mix = strings(height(Tb), 1);
            if ismember("Mixture", vn), mix = strtrim(string(Tb.Mixture)); mix(ismissing(mix)) = ""; end
            bad = ~ismember(b, ["psig", "psia"]) & mix ~= "";
            if any(bad)
                error('PressureBasis:BadFlag', 'Sheet "%s", data row(s) %s: %s must be "psig" or "psia".', ...
                    sheet, strjoin(string(find(bad).'), ", "), entities.PressureBasis.FlagColumn);
            end
            off(b == "psig") = entities.PressureBasis.Atm_psi;
        end

        function Tb = apply(Tb, mode, sheet)
            % Converts the pressure columns of Tb to the requested basis and adds a
            % column P_offset_psi with the offset applied to each row.
            arguments
                Tb    table
                mode  (1,1) string {mustBeMember(mode, ["absolute", "recorded"])}
                sheet (1,1) string = "(table)"
            end
            off = entities.PressureBasis.offset(Tb, mode, sheet);
            vn = string(Tb.Properties.VariableNames);
            for c = entities.PressureBasis.Columns(ismember(entities.PressureBasis.Columns, vn))
                v = Tb.(c);
                if ~isnumeric(v), v = str2double(string(v)); end
                Tb.(c) = double(v) + off;
            end
            Tb.P_offset_psi = off;
        end

        function Tb = readSheet(file, sheet, mode)
            % readtable (names preserved) followed by apply.
            arguments
                file  (1,1) string
                sheet (1,1) string
                mode  (1,1) string {mustBeMember(mode, ["absolute", "recorded"])}
            end
            Tb = readtable(file, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
            Tb = entities.PressureBasis.apply(Tb, mode, sheet);
        end
    end
end
