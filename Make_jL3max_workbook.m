%% MAKE_JL3MAX_WORKBOOK  Workbook for the upper-BIP calibration jL3-max (variant "psiamax")
%
% jL3-max = jL3 with each C1 pair BIP at the upper edge of its measurement-based
% window and the four group TIPs refitted at those BIPs (corner "HHH" of
% Uncertainty_BIP_TIP.m). k, p, lambda, the nC5-nC10 BIP, the late B1 radius and
% everything else are those of jL3.
%
% Reads   results/uncertainty_bip_tip_psia/uncertainty_bip_tip.xlsx
%           sheet Windows   : Window_hi at CompRel_pct = cfg.compRel_pct
%           sheet GroupTIPs : TIP at Level 3 (high edge), same CompRel_pct
%         config/MixtureData - Recalibrated psia.xlsx   (jL3 workbook)
% Writes  config/MixtureData - Recalibrated psia max.xlsx
%           sheets BIP and TIP_Groups replaced; a row "jL3max" added to Calibration_Info
%
% Then run MLDO_check.m, Temperature_test.m and NetStress_test.m (and PVT_verification.m
% if the GEM export is wanted) with cfg.variant = "psiamax". The preset is "jL3", which
% reads TIP_Groups from this workbook, so its rock-offset correlation is the HHH one.
%
% Needs (local): the two files above (run Uncertainty_BIP_TIP.m with CompRel [0.05 0.10]
% first).

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

cfg.compRel_pct = 10;                                          % composition scenario of jL3-max
cfg.source = fullfile(scriptDir, "config", "MixtureData - Recalibrated psia.xlsx");
cfg.target = fullfile(scriptDir, "config", "MixtureData - Recalibrated psia max.xlsx");
cfg.unc    = fullfile(scriptDir, "results", "uncertainty_bip_tip_psia", "uncertainty_bip_tip.xlsx");
cfg.pairs  = ["C1-nC5", "C1-nC8", "C1-nC10"];

if ~isfile(cfg.source), error('jl3max:NoWorkbook', 'Not found: %s', cfg.source); end
if ~isfile(cfg.unc),    error('jl3max:NoUncertainty', 'Not found: %s (run Uncertainty_BIP_TIP.m)', cfg.unc); end

% Upper window edges and group TIPs at the high level
W = readtable(cfg.unc, 'Sheet', 'Windows', 'TextType', 'string');
G = readtable(cfg.unc, 'Sheet', 'GroupTIPs', 'TextType', 'string');
bipMax = nan(1, numel(cfg.pairs));
for n = 1:numel(cfg.pairs)
    r = find(W.Pair == cfg.pairs(n) & W.CompRel_pct == cfg.compRel_pct, 1);
    if isempty(r), error('jl3max:NoWindow', 'Windows has no row %s at %g %%.', cfg.pairs(n), cfg.compRel_pct); end
    if W.Window_lo(r) > W.Window_hi(r), error('jl3max:Empty', 'Empty window for %s.', cfg.pairs(n)); end
    bipMax(n) = W.Window_hi(r);
end
groups = thermo.RockOffsetTIP.Groups;
tipMax = nan(1, numel(groups));
for g = 1:numel(groups)
    r = find(strtrim(G.Group) == groups(g) & G.CompRel_pct == cfg.compRel_pct & G.Level == 3, 1);
    if isempty(r), error('jl3max:NoTIP', 'GroupTIPs has no level-3 row for %s.', groups(g)); end
    tipMax(g) = G.TIP(r);
end

copyfile(cfg.source, cfg.target, 'f');

% Sheet BIP (symmetric)
C = readcell(cfg.target, 'Sheet', 'BIP');
miss = cellfun(@(e) isa(e, 'missing'), C);  C(miss) = {''};
hdr = string(C(1, 2:end));  lab = string(C(2:end, 1));
for n = 1:numel(cfg.pairs)
    pr = split(cfg.pairs(n), "-");
    ra = find(lab == pr(1), 1);  rb = find(lab == pr(2), 1);  ca = find(hdr == pr(1), 1);  cb = find(hdr == pr(2), 1);
    C{1 + ra, 1 + cb} = bipMax(n);  C{1 + rb, 1 + ca} = bipMax(n);
end
writecell(C, cfg.target, 'Sheet', 'BIP', 'Range', 'A1');

% Sheet TIP_Groups (same layout as Recalibrate_psia.m writes)
TG = readtable(cfg.target, 'Sheet', 'TIP_Groups', 'TextType', 'string', 'VariableNamingRule', 'preserve');
for g = 1:numel(groups)
    r = find(strtrim(TG.Group) == groups(g), 1);
    TG.TIP(r) = tipMax(g);
end
if ismember('Objective', TG.Properties.VariableNames)
    TG.Objective(:) = "jL3-max: refitted at the " + cfg.compRel_pct + " % upper BIP edge (Uncertainty_BIP_TIP.m)";
end
writetable(TG, cfg.target, 'Sheet', 'TIP_Groups');

% Calibration_Info: note (late B1 radius kept from jL3)
CI = readtable(cfg.target, 'Sheet', 'Calibration_Info', 'TextType', 'string');
CI = [CI; {"jL3max", sprintf('upper BIP edge at %g %%: C1-nC5 %.6f, C1-nC8 %.6f, C1-nC10 %.6f; late B1 radius of jL3', ...
            cfg.compRel_pct, bipMax)}];
writetable(CI, cfg.target, 'Sheet', 'Calibration_Info');

% Check: the jL3 preset on the new workbook reproduces the HHH correlation
tip = thermo.RockOffsetTIP.frozen("jL3", cfg.target);
fprintf('jL3-max written: %s\n', cfg.target);
fprintf('  BIP  C1-nC5 %.6f | C1-nC8 %.6f | C1-nC10 %.6f\n', bipMax);
fprintf('  TIPs %s\n', strjoin(compose("%.4f", tipMax), " / "));
fprintf('  correlation alpha_EF2 %.5f | alpha_B1 %.5f | B %.5f\n', tip.alphaEF2, tip.alphaB1, tip.B);
