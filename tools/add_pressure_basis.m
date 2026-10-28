function report = add_pressure_basis(file, opts)
% ADD_PRESSURE_BASIS  Write the P_basis column ("psig" | "psia") into the pressure
% sheets of a mixture workbook (in place). The recorded values are not changed;
% entities.PressureBasis converts them when a script loads with "absolute".
%
%   report = add_pressure_basis("config/MixtureData - Recalibrated psia.xlsx")
%
% Rule (laboratory answer of 28 Oct 2026; Doryani thesis tables = early campaign):
%   sheet with a Campaign column : early -> psig (LabView VI / Druck display unit)
%                                  late  -> psia (Quartzdyne with Q-Com, from ~2017)
%   sheet without Campaign       : Source naming the Doryani thesis -> psig
%   anything else                : error (state the basis by hand)
% An existing P_basis column is kept unless opts.Overwrite is true.
%
% report: table Sheet | Rows | psig | psia | Action

arguments
    file (1,1) string
    opts.Sheets (1,:) string = ["Mix_Props", "Bulk_PVT", "Mix_Props_T", "Mix_Props_NS"]
    opts.Overwrite (1,1) logical = false
end
if ~isfile(file), error('pbasis:NoFile', 'Not found: %s', file); end
flag = entities.PressureBasis.FlagColumn;
have = string(sheetnames(file));
rows = {};
for s = opts.Sheets
    if ~ismember(s, have), rows(end+1, :) = {s, 0, 0, 0, "sheet absent"}; continue; end %#ok<AGROW>
    C = readcell(file, 'Sheet', s);
    miss = cellfun(@(e) isa(e, 'missing'), C);
    C(miss) = {''};
    hdr = string(C(1, :));
    iMix = find(hdr == "Mixture", 1);
    if isempty(iMix), error('pbasis:NoMixture', 'Sheet %s has no Mixture column.', s); end
    data = find(arrayfun(@(r) strlength(strtrim(string(C{r, iMix}))) > 0, 2:size(C, 1))) + 1;
    iFlag = find(hdr == flag, 1);
    if ~isempty(iFlag) && ~opts.Overwrite
        b = lower(strtrim(string(C(data, iFlag))));
        rows(end+1, :) = {s, numel(data), nnz(b == "psig"), nnz(b == "psia"), "kept existing column"}; %#ok<AGROW>
        continue;
    end
    iCamp = find(lower(hdr) == "campaign", 1);
    iSrc  = find(lower(hdr) == "source", 1);
    b = strings(numel(data), 1);
    for n = 1:numel(data)
        r = data(n);
        if ~isempty(iCamp) && any(lower(strtrim(string(C{r, iCamp}))) == ["early", "late"])
            if lower(strtrim(string(C{r, iCamp}))) == "early", b(n) = "psig"; else, b(n) = "psia"; end
        elseif ~isempty(iSrc) && contains(string(C{r, iSrc}), "Doryani", 'IgnoreCase', true)
            b(n) = "psig";
        else
            error('pbasis:Unknown', ['Sheet %s, row %d: no early/late Campaign and no Doryani Source; ' ...
                'enter %s by hand.'], s, r, flag);
        end
    end
    if isempty(iFlag), iFlag = size(C, 2) + 1; end
    col = cell(size(C, 1), 1);  col(:) = {''};
    col{1} = char(flag);
    col(data) = cellstr(b);
    writecell(col, file, 'Sheet', s, 'Range', localColumn(iFlag) + "1");
    rows(end+1, :) = {s, numel(data), nnz(b == "psig"), nnz(b == "psia"), "written (column " + localColumn(iFlag) + ")"}; %#ok<AGROW>
end
report = cell2table(rows, 'VariableNames', {'Sheet', 'Rows', 'psig', 'psia', 'Action'});
end

function L = localColumn(n)
% Excel column letters of column number n.
L = "";
while n > 0
    m = mod(n - 1, 26);
    L = char('A' + m) + L;
    n = floor((n - 1) / 26);
end
end
