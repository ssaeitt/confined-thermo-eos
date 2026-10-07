%% BACKCALC_RP_B1_LATE  Pore-radius range of the late B1 core from C1-nC8 [0.9, 0.1]
%
% Why: the late-core B1 C1-nC8 10% point (bulk 4126 psi, confined 4342-4361 psi) cannot be
% reconciled by adjusting the composition, so the pore radius is back-calculated from the
% EXPERIMENTAL CONFINED Pdew range instead.
%
% Model (frozen, FWI-Only; cfg.preset selects joint-Linearized or pareto-exact-w1):
%   a_eff,ij = sqrt(ai aj)(1-kij) - sqrt(ci cj)(1-TIP_ij)
%   TIP_ij   = -exp( alpha_rock + B * x_ij )          (rock-offset, size-independent)
%   primary      LL-inv_wbar   x_ij = 1 / wbar_ij,  wbar_ij = (omega_i + omega_j)/2
%   sensitivity  LL-inv_Tstar  x_ij = sqrt(eps_i eps_j) / (kB T)  (= 1/T*)
% The TIP does not depend on r_p, so r_p acts only through c_i and the shift must fall
% monotonically as r_p rises (checked below).
%
% Targets (both reported; the first is the primary):
%   absolute        model P_dew,conf(r_p) inside the experimental confined range
%   shift_modelBulk model shift (conf - MODEL bulk) inside the experimental shift range
% The r_p range is [r_p(P = upper), r_p(P = lower)]; a higher Pdew needs a smaller pore.
%
% Needs: config/MixtureData - Complete Mixtures.xlsx, config/tip_correlation_rockoffset.xlsx

clear; clc; close all;

%% Preflight ----------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
if isfolder(fullfile(scriptDir, "src")), addpath(fullfile(scriptDir, "src")); end
for pk = ["entities", "thermo", "solvers"]
    if isempty(meta.package.fromName(pk))
        cand = dir(fullfile(scriptDir, "**", "+" + pk));
        cand = unique(string({cand.folder}));
        cand = cand(endsWith(cand, filesep + "+" + pk));
        if isempty(cand)
            error('bc:PackageNotFound', 'Folder "+%s" not found under %s.', pk, scriptDir);
        end
        addpath(fileparts(cand(1)));
    end
end

%% 0. Configuration -------------------------------------------------------------
cfg.configDir   = fullfile(scriptDir, "config");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.mixtureFile = fullfile(cfg.configDir, "MixtureData - Complete Mixtures.xlsx");
cfg.tipFile     = fullfile(cfg.configDir, "tip_correlation_rockoffset.xlsx");
cfg.outFile     = "backcalc_rp_B1_late";
cfg.T           = 293.15;
cfg.mode        = "FWIOnly";
cfg.kB          = 1.380649e-23;                        % [J/K]
cfg.poreRadiusNominal = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});

% Parameter preset. jL = joint-Linearized (primary, inv_wbar / inv_Tstar).
% w1 = pareto-exact-w1 with the Tstar forms that passed the criteria (sensitivity); its results are
% saved under a "_w1" suffix so they never overwrite the jL file.
cfg.preset = "w1";
presets.jL = struct('combination', "joint-linearized", 'k', 112.42634826677, 'p', 1.3558179290678, ...
    'lambda', 2.0916235090532, 'forms', ["LL-inv_wbar", "LL-inv_Tstar"], 'suffix', "");
presets.w1 = struct('combination', "pareto-exact-w1", 'k', 166.410385934507, 'p', 1.5488421885245, ...
    'lambda', 4.47908205631447, 'forms', ["LL-Tstar", "LL-ln_Tstar"], 'suffix', "_w1");
ps = presets.(cfg.preset);
cfg.fwi = struct('k', ps.k, 'p', ps.p, 'lambda', ps.lambda);
cfg.outFile = cfg.outFile + ps.suffix;

% TIP correlation rows to use (sheet "Correlation")
cfg.tip.combination = ps.combination;
cfg.tip.mode        = "FWIOnly";
cfg.models = struct('role', {"primary", "sensitivity"}, 'form', {ps.forms(1), ps.forms(2)});

% The case that sets r_p: late B1 C1-nC8 10 % nC8 (identified by composition)
cfg.rp.case       = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.rp.bracket_nm = [4, 12];        % search range
cfg.rp.step_nm    = 0.1;            % scan step (monotonicity check and root bracketing)
cfg.rp.band_nm    = [6, 11.25];     % reporting band (metadata)
cfg.rp.tolX_nm    = 1e-3;
cfg.maxShift_psi  = Inf;            % lost-branch guard off: small pores may exceed 500 psi
cfg.gateTol       = 1e-6;           % TIP regression gate vs the Pair_TIPs sheet

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
psi = entities.DewPointDataset.psi2Pa;

%% 1. Data and engine ------------------------------------------------------------
dp  = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
idx = dp.findCases(cfg.rp.case.Mixture, cfg.rp.case.Rock, cfg.rp.case.z);
if numel(idx) ~= 1
    error('bc:CaseNotUnique', '%d rows match %s / %s z = [%s]; expected exactly one.', numel(idx), ...
        cfg.rp.case.Mixture, cfg.rp.case.Rock, num2str(cfg.rp.case.z.'));
end
if ~dp.HasConfined(idx)
    error('bc:NoConfinedData', 'The selected row has no confined Pdew range.');
end
c = dp.getCase(idx);
exp_ = struct('bulk', dp.BulkPdew_psi(idx), 'lo', dp.ConfLower_psi(idx), 'up', dp.ConfUpper_psi(idx));
exp_.mid = 0.5 * (exp_.lo + exp_.up);
exp_.shLo = exp_.lo - exp_.bulk;  exp_.shUp = exp_.up - exp_.bulk;
fprintf('Case: %s / %s, z = [%s], T = %.2f K\n', c.Mixture, c.Rock, num2str(c.z.'), c.T);
fprintf('Experiment: bulk %.0f psi | confined [%.0f, %.0f] psi | shift [%.0f, %.0f] psi\n\n', ...
    exp_.bulk, exp_.lo, exp_.up, exp_.shLo, exp_.shUp);

runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadiusNominal, 'maxShift_psi', cfg.maxShift_psi);
eng = runner.engine(c.Mixture, c.Rock);
runner.setFWI(cfg.fwi.k, cfg.fwi.p, cfg.fwi.lambda);
fluid = eng.eos.Fluid;

[Pb_Pa, ~, okBulk] = runner.bulk(c);
if ~okBulk, error('bc:BulkFailed', 'Model bulk dew point did not converge.'); end
Pb = Pb_Pa / psi;
fprintf('Model bulk Pdew (literature BIPs, charged composition): %.2f psi (exp %.0f, %+.2f %%)\n\n', ...
    Pb, exp_.bulk, 100 * (Pb / exp_.bulk - 1));

% Targets on the two bases
tg.absolute        = struct('lo', exp_.lo,      'up', exp_.up,      'mid', exp_.mid);
tg.shift_modelBulk = struct('lo', Pb + exp_.shLo, 'up', Pb + exp_.shUp, 'mid', Pb + 0.5 * (exp_.shLo + exp_.shUp));
bases = ["absolute", "shift_modelBulk"];

%% 2. TIP from the rock-offset correlation (+ regression gate) -------------------
tipTab  = readtable(cfg.tipFile, 'Sheet', 'Correlation', 'VariableNamingRule', 'preserve');
pairTab = readtable(cfg.tipFile, 'Sheet', 'Pair_TIPs',   'VariableNamingRule', 'preserve');

nm = string(fluid.ComponentNames);
pairName = nm(1) + "-" + nm(2);
tipSet = cell(1, numel(cfg.models));
for m = 1:numel(cfg.models)
    co = tipCoefficients(tipTab, cfg.tip.combination, cfg.tip.mode, cfg.models(m).form, c.Rock);
    X  = pairDescriptor(fluid, co.descriptor, cfg.T, cfg.kB);
    K  = -exp(co.alpha + co.B * X);
    K(1:numel(nm)+1:end) = 0;
    K  = (K + K.') / 2;

    % Gate: the pair TIP must reproduce the Pair_TIPs sheet for this rock
    sel = string(pairTab.Combination) == cfg.tip.combination & string(pairTab.Mode) == cfg.tip.mode & ...
          string(pairTab.Form) == cfg.models(m).form & string(pairTab.Context) == "binary" & ...
          string(pairTab.Rock) == c.Rock & string(pairTab.Pair) == pairName;
    if nnz(sel) ~= 1
        error('bc:GateRow', 'Pair_TIPs row not found/unique for %s %s %s.', cfg.models(m).form, c.Rock, pairName);
    end
    ref = pairTab.TIP(sel);
    if abs(K(1, 2) - ref) > cfg.gateTol * max(1, abs(ref))
        error('bc:TIPGate', '%s: TIP(%s, %s) = %.8f vs sheet %.8f.', cfg.models(m).form, c.Rock, pairName, K(1, 2), ref);
    end
    tipSet{m} = K;
    fprintf('%-12s (%s): alpha_%s = %+.6f, B = %+.6f, x = %.6f -> TIP(%s) = %+.6f  [gate vs sheet OK]\n', ...
        cfg.models(m).form, cfg.models(m).role, c.Rock, co.alpha, co.B, X(1, 2), pairName, K(1, 2));
end
if numel(tipSet) == 2 && abs(tipSet{1}(1, 2) - tipSet{2}(1, 2)) < 1e-6
    fprintf(['Note: both forms give the same B1 %s TIP (alpha_B1 is fitted through the single B1 binary), so ' ...
             'the r_p back-calculation is identical for primary and sensitivity; they differ only for other pairs.\n'], pairName);
end
fprintf('\n');

%% 3. r_p scan and root finding ---------------------------------------------------
grid = cfg.rp.bracket_nm(1):cfg.rp.step_nm:cfg.rp.bracket_nm(2);
if ~any(abs(grid - 11.25) < 1e-9), grid = sort([grid, 11.25]); end
scanRows = {};  resRows = {};  derived = struct();

for m = 1:numel(cfg.models)
    form = cfg.models(m).form;  K = tipSet{m};
    runner.clearWarm();
    fun = @(rnm) pConf(rnm, runner, c, K, form, cfg.mode, psi);

    P = nan(size(grid));  route = strings(size(grid));
    for g = 1:numel(grid)
        [P(g), route(g)] = fun(grid(g));
        scanRows(end+1, :) = {form, grid(g), P(g), P(g) - Pb, route(g)}; %#ok<AGROW>
    end
    fin = isfinite(P);
    monotone = nnz(fin) >= 3 && all(diff(P(fin)) < 0);
    fprintf('--- %s (%s): %d/%d scan points solved | P_dew falls monotonically with r_p: %s ---\n', ...
        form, cfg.models(m).role, nnz(fin), numel(grid), mat2str(monotone));
    if ~monotone
        warning('bc:NotMonotone', '%s: P_dew(r_p) is not strictly decreasing on the scan; check the scan table.', form);
    end

    [~, jNom] = min(abs(grid - 11.25));
    fprintf('    at nominal 11.25 nm: P_dew = %.2f psi (shift %+.2f psi vs model bulk)\n', P(jNom), P(jNom) - Pb);

    for b = bases
        t = tg.(b);
        rUp  = pickRoot(fun, grid, P, t.up,  cfg.rp.tolX_nm, cfg.rp.band_nm);   % higher P -> smaller r_p
        rMid = pickRoot(fun, grid, P, t.mid, cfg.rp.tolX_nm, cfg.rp.band_nm);
        rLo  = pickRoot(fun, grid, P, t.lo,  cfg.rp.tolX_nm, cfg.rp.band_nm);
        closure = NaN;
        if isfinite(rMid.value), closure = fun(rMid.value) - t.mid; end
        inBand = isfinite(rMid.value) && rMid.value >= cfg.rp.band_nm(1) && rMid.value <= cfg.rp.band_nm(2);
        resRows(end+1, :) = {form, cfg.models(m).role, b, t.lo, t.up, rUp.value, rMid.value, rLo.value, ...
            0.5 * (rLo.value - rUp.value), max(rUp.n, max(rMid.n, rLo.n)), inBand, closure, monotone}; %#ok<AGROW>
        fprintf('    %-16s target P [%.1f, %.1f] psi -> r_p = %.3f nm  (range %.3f to %.3f nm) | roots: %d | in band: %s | closure %+.3f psi\n', ...
            b, t.lo, t.up, rMid.value, rUp.value, rLo.value, max(rUp.n, max(rMid.n, rLo.n)), mat2str(inBand), closure);
        if b == "absolute"
            derived.(matlab.lang.makeValidName(form)) = struct('rp_mid_m', rMid.value * 1e-9, ...
                'rp_lo_m', rUp.value * 1e-9, 'rp_hi_m', rLo.value * 1e-9);
        end
    end
    fprintf('\n');
end

%% 4. Tables, files ----------------------------------------------------------------
scanTable = textToString(cell2table(scanRows, 'VariableNames', {'Form', 'rp_nm', 'Pdew_conf_psi', ...
    'Shift_modelBulk_psi', 'Route'}));
resTable = textToString(cell2table(resRows, 'VariableNames', {'Form', 'Role', 'Basis', 'Target_lo_psi', ...
    'Target_up_psi', 'rp_at_Pup_nm', 'rp_mid_nm', 'rp_at_Plo_nm', 'HalfRange_nm', 'NumRoots', 'InBand', ...
    'Closure_psi', 'MonotoneScan'}));
fprintf('=== Summary ===\n'); disp(resTable);

info = struct('case', c.Mixture + " " + c.Rock + " z=[" + strjoin(string(c.z.'), ",") + "]", ...
    'expBulk_psi', exp_.bulk, 'expConf_psi', [exp_.lo, exp_.up], 'modelBulk_psi', Pb, ...
    'tipFile', cfg.tipFile, 'combination', cfg.tip.combination, 'mode', cfg.tip.mode, 'fwi', cfg.fwi);
out = fullfile(cfg.resultsDir, cfg.outFile);
save(out + ".mat", 'cfg', 'info', 'resTable', 'scanTable', 'derived');
xls = out + ".xlsx";
if isfile(xls), delete(xls); end
writetable(resTable,  xls, 'Sheet', 'rp_result');
writetable(scanTable, xls, 'Sheet', 'rp_scan');
fprintf('Written: %s\n', xls);
fprintf('Derived r_p (absolute basis) saved in %s (variable "derived", metres).\n', out + ".mat");

%% ============================== Local functions ==============================

function co = tipCoefficients(tab, comb, mode, form, rock)
% Coefficients of ln(-TIP) = alpha_rock + B * x from the "Correlation" sheet.
sel = string(tab.Combination) == comb & string(tab.Mode) == mode & string(tab.Form) == form;
if nnz(sel) ~= 1
    error('bc:CorrelationRow', '%d rows match %s / %s / %s in the Correlation sheet.', nnz(sel), comb, mode, form);
end
vn = "alpha_" + rock;
if ~ismember(vn, string(tab.Properties.VariableNames))
    error('bc:NoRockOffset', 'Column %s not found: no offset for rock %s.', vn, rock);
end
co = struct('alpha', tab.(vn)(sel), 'B', tab.B(sel), 'descriptor', string(tab.Descriptor(sel)));
end

function X = pairDescriptor(fluid, name, T, kB)
% Dimensionless pair descriptors (match the x column of the Pair_TIPs sheet).
w  = fluid.omega(:);
e  = fluid.LJ_Energy(:);                      % [J]
MW = fluid.MW(:);
Ts = (kB * T) ./ sqrt(e * e.');               % T* = kB T / sqrt(eps_i eps_j)
switch name
    case "wbar",      X = 0.5 * (w + w.');
    case "inv_wbar",  X = 2 ./ (w + w.');
    case "Tstar",     X = Ts;
    case "inv_Tstar", X = 1 ./ Ts;
    case "ln_Tstar",  X = log(Ts);
    case "inv_Mgeo",  X = 1 ./ sqrt(MW * MW.');
    otherwise
        error('bc:UnknownDescriptor', 'Descriptor "%s" is not implemented.', name);
end
end

function [Ppsi, route] = pConf(rnm, runner, c, K, form, mode, psi)
% Confined Pdew [psi] at r_p = rnm [nm] with the given TIP matrix; NaN if no accepted root.
r = runner.solve(c, mode, K, 'r', rnm * 1e-9, 'Key', form + "|" + sprintf('%.5f', rnm));
Ppsi = NaN;  route = "FAIL: " + r.reason;
if r.ok, Ppsi = r.P / psi;  route = string(r.route); end
end

function out = pickRoot(fun, grid, P, target, tol, band)
% All crossings of P_dew(r_p) = target, refined with fzero; choose the one inside the band
% (or the nearest to it). value = NaN if the target is not bracketed on the scan.
roots = [];
g = P - target;
for j = 1:numel(grid) - 1
    if isfinite(g(j)) && isfinite(g(j+1)) && g(j) * g(j+1) < 0
        a = grid(j);  b = grid(j+1);
        try
            x = fzero(@(x) fun(x) - target, [a b], optimset('TolX', tol));
        catch
            x = a - g(j) * (b - a) / (g(j+1) - g(j));
        end
        roots(end+1) = x; %#ok<AGROW>
    end
end
out = struct('value', NaN, 'n', numel(roots), 'all', roots);
if isempty(roots), return; end
d = arrayfun(@(x) max([band(1) - x, 0, x - band(2)]), roots);
[~, j] = min(d);
out.value = roots(j);
end

function t = textToString(t)
for v = string(t.Properties.VariableNames)
    col = t.(v);
    if iscell(col) && ~isempty(col) && all(cellfun(@(e) isstring(e) || ischar(e), col))
        t.(v) = string(col);
    end
end
end