%% NETSTRESS_TEST  Blind test of the pore-radius dependence with Doryani's net-stress data
%
% Data: sheet Mix_Props_NS. Four fluids on the EF2 core (thesis IDs F4, F5, F7, F10;
% labels C5(5), C10(8), C8(8), C5C10(4,4)) at net stress 500, 1000 and 2000 psi
% (Doryani thesis Sec. 4.2.3). Raising the net stress compresses the core, so the same
% fluid is measured on the same rock at three mean pore radii.
%
% Pore radius per net stress (two sets, both reported):
%   "thesis" : 54 / 45 / 38 nm (Doryani Sec. 5.5.3, from Nazari Moghaddam & Jamiolahmady
%              2016, pore-size distribution scaled uniformly with stress). Primary set.
%   "sqrtk"  : r = 54 nm * sqrt(k_ref / k_ref,500), per fluid (Kozeny-type scaling at
%              constant porosity). Sensitivity check only.
%
% Model: frozen FWI-Only preset (jL, or jL2 for the recalibration). Nothing is refitted:
% the rock-offset TIP alpha_EF2 belongs to the rock chemistry, which stress does not
% change, so only r enters through the wall term s(sigma/r). The 500 psi rows are the
% calibration points; 1000 and 2000 psi are blind.
%
% Reports, per row and radius set: model bulk and confined Pdew, shift vs model bulk,
% measured shift range, in-range flags; per row, the pore radius that reproduces the
% measured lower / mid / upper confined Pdew (back-calculated radius); per fluid, a
% dense radius sweep of the model shift.

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
% --- Model variant -----------------------------------------------------------------
%   "frozen" : current calibration (workbook "MixtureData - Complete Mixtures.xlsx", preset jL)
%   "recal"  : C1-nC5 recalibration (workbook "MixtureData - Recalibrated.xlsx", preset jL2)
cfg.variant = "recal";
variants = struct( ...
    'frozen', struct('file', "MixtureData - Complete Mixtures.xlsx", 'preset', "jL"), ...
    'recal',  struct('file', "MixtureData - Recalibrated.xlsx",     'preset', "jL2"));
vsel = variants.(char(cfg.variant));
cfg.mixtureFile = fullfile(scriptDir, "config", vsel.file);
cfg.preset      = vsel.preset;
cfg.sheet       = "Mix_Props_NS";
cfg.outDir      = fullfile(scriptDir, "results", "netstress_test_" + cfg.variant);
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
cfg.radiusSets  = ["thesis", "sqrtk"];
cfg.rRef_nm     = 54;                          % EF2 mean pore radius at 500 psi net stress
cfg.rSweep_nm   = [15 20 25 30 34 38 42 45 50 54 60 70 80 100 130];
cfg.backCalc    = true;
cfg.rBracket_nm = [2, 1000];                   % search range for the back-calculated radius

if ~isfile(cfg.mixtureFile)
    error('nstest:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.mixtureFile);
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'Sheet', cfg.sheet);
mt = readtable(cfg.mixtureFile, 'Sheet', cfg.sheet, 'VariableNamingRule', 'preserve');
fid   = string(mt.Fluid_ID);
label = string(mt.Label);
ns    = mt.Net_stress_psi;
kref  = mt.k_ref_mD;
if any(~isfinite(dp.PoreRadius_nm))
    error('nstest:NoRadius', 'Sheet %s needs an rp_nm value on every row.', cfg.sheet);
end
fprintf('%d rows in %s; preset %s\n\n', dp.NumCases, cfg.sheet, cfg.preset);

tip = thermo.RockOffsetTIP.frozen(cfg.preset, cfg.mixtureFile);
fprintf('%s\n\n', tip.summary());
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius, 'maxShift_psi', Inf);
runner.setFWI(tip.k, tip.p, tip.lambda);

% Radius of every row under each radius set [nm]
rSet = struct();
rSet.thesis = dp.PoreRadius_nm;
rSet.sqrtk  = nan(dp.NumCases, 1);
for f = unique(fid).'
    sel = fid == f;
    k500 = kref(sel & ns == 500);
    rSet.sqrtk(sel) = cfg.rRef_nm * sqrt(kref(sel) / k500(1));
end

%% 1. Model at the measured net stresses -------------------------------------------
rows = {};
for i = 1:dp.NumCases
    c   = dp.getCase(i);
    eng = runner.engine(c.Mixture, c.Rock);
    K   = tip.matrix(eng.eos.Fluid, c.Rock, c.T);
    lo  = dp.ConfLower_psi(i);  up = dp.ConfUpper_psi(i);  Pbe = dp.BulkPdew_psi(i);
    role = "blind";  if ns(i) == 500, role = "calibration"; end

    for rs = cfg.radiusSets
        r_nm = rSet.(char(rs))(i);
        res  = runner.solve(c, "FWIOnly", K, 'r', r_nm * 1e-9, 'Key', "ns|" + fid(i));
        [Pc, Pbm, shM] = deal(NaN);
        if res.ok
            Pc = res.P / psi;  Pbm = res.PdewBulkModel / psi;  shM = Pc - Pbm;
        end
        rows(end+1, :) = {fid(i), label(i), ns(i), role, rs, r_nm, Pbe, Pbm, lo, up, Pc, ...
            100 * (Pc / (0.5 * (lo + up)) - 1), Pc >= lo && Pc <= up, ...
            lo - Pbe, up - Pbe, shM, shM >= lo - Pbe && shM <= up - Pbe, string(res.reason)}; %#ok<AGROW>
        fprintf('  %-11s %4d psi  %-6s r = %5.1f nm | conf %7.1f (exp %4.0f-%4.0f) | shift %6.1f (exp %3.0f-%3.0f)\n', ...
                label(i), ns(i), rs, r_nm, Pc, lo, up, shM, lo - Pbe, up - Pbe);
    end
end
R = cell2table(rows, 'VariableNames', {'Fluid_ID', 'Label', 'Net_stress_psi', 'Role', 'RadiusSet', ...
    'r_nm', 'Bulk_exp_psi', 'Bulk_model_psi', 'Conf_exp_lo_psi', 'Conf_exp_up_psi', 'Conf_model_psi', ...
    'Conf_err_vs_mid_pct', 'Conf_in_range', 'Shift_exp_lo_psi', 'Shift_exp_up_psi', ...
    'Shift_model_psi', 'Shift_in_range', 'Reason'});
R = convertvars(R, {'Fluid_ID', 'Label', 'Role', 'RadiusSet', 'Reason'}, 'string');
disp(R);
writetable(R, fullfile(cfg.outDir, "netstress_test.csv"));

%% 2. Back-calculated radius: r that reproduces the measured confined Pdew -----------
if cfg.backCalc
    bc = {};
    for i = 1:dp.NumCases
        c   = dp.getCase(i);
        eng = runner.engine(c.Mixture, c.Rock);
        K   = tip.matrix(eng.eos.Fluid, c.Rock, c.T);
        tgt = [dp.ConfUpper_psi(i), 0.5 * (dp.ConfLower_psi(i) + dp.ConfUpper_psi(i)), dp.ConfLower_psi(i)];
        rb  = nan(1, 3);
        for n = 1:3
            rb(n) = localBackCalc(runner, c, K, tgt(n) * psi, cfg.rBracket_nm, "bc|" + fid(i));
        end
        bc(end+1, :) = {fid(i), label(i), ns(i), dp.PoreRadius_nm(i), rSet.sqrtk(i), rb(1), rb(2), rb(3)}; %#ok<AGROW>
        fprintf('  %-11s %4d psi: back-calculated r = %6.1f nm (range %6.1f - %6.1f nm); thesis %4.1f, sqrt(k) %4.1f\n', ...
                label(i), ns(i), rb(2), rb(1), rb(3), dp.PoreRadius_nm(i), rSet.sqrtk(i));
    end
    B = cell2table(bc, 'VariableNames', {'Fluid_ID', 'Label', 'Net_stress_psi', 'r_thesis_nm', ...
        'r_sqrtk_nm', 'r_back_at_upper_nm', 'r_back_at_mid_nm', 'r_back_at_lower_nm'});
    B = convertvars(B, {'Fluid_ID', 'Label'}, 'string');
    disp(B);
    writetable(B, fullfile(cfg.outDir, "netstress_backcalc_radius.csv"));
end

%% 3. Dense radius sweep per fluid ---------------------------------------------------------
sw = {};
[~, first] = unique(fid, 'stable');
for i = first.'
    c   = dp.getCase(i);
    eng = runner.engine(c.Mixture, c.Rock);
    K   = tip.matrix(eng.eos.Fluid, c.Rock, c.T);
    for r_nm = sort(cfg.rSweep_nm, 'descend')          % large -> small r: warm start follows the shift
        res = runner.solve(c, "FWIOnly", K, 'r', r_nm * 1e-9, 'Key', "sweep|" + fid(i));
        sh = NaN;  if res.ok, sh = (res.P - res.PdewBulkModel) / psi; end
        sw(end+1, :) = {fid(i), label(i), r_nm, sh, string(res.reason)}; %#ok<AGROW>
    end
end
S = cell2table(sw, 'VariableNames', {'Fluid_ID', 'Label', 'r_nm', 'Shift_model_psi', 'Reason'});
S = convertvars(S, {'Fluid_ID', 'Label', 'Reason'}, 'string');
writetable(S, fullfile(cfg.outDir, "netstress_radius_sweep.csv"));

%% 4. Summary: change in shift from 500 to 2000 psi -----------------------------------------
fprintf('\nChange in shift, 500 -> 2000 psi net stress (thesis radii 54 -> 38 nm):\n');
for i = first.'
    a = R(R.Fluid_ID == fid(i) & R.RadiusSet == "thesis" & R.Net_stress_psi == 500, :);
    b = R(R.Fluid_ID == fid(i) & R.RadiusSet == "thesis" & R.Net_stress_psi == 2000, :);
    dExp = 0.5 * (b.Shift_exp_lo_psi + b.Shift_exp_up_psi) - 0.5 * (a.Shift_exp_lo_psi + a.Shift_exp_up_psi);
    fprintf('  %-11s measured %+6.1f psi | model %+6.1f psi\n', label(i), dExp, b.Shift_model_psi - a.Shift_model_psi);
end

%% 5. Figure: shift vs pore radius, one tile per fluid ----------------------------------------
figure('Name', 'Net-stress test', 'Position', [60 60 1150 780]);
tiledlayout(2, 2, 'TileSpacing', 'compact');
% Plotted over the measured radius window only (the sweep CSV keeps the full range); the
% sqrt(k) radius set is reported in the CSVs but not plotted (rejected by the test).
cfg.plotR_nm = [30, 80];
for i = first.'
    nexttile; hold on; box on; grid on;
    s = S(S.Fluid_ID == fid(i) & S.r_nm >= cfg.plotR_nm(1) & S.r_nm <= cfg.plotR_nm(2), :);
    s = sortrows(s, 'r_nm');
    sel = find(fid == fid(i));
    lo = dp.ConfLower_psi(sel) - dp.BulkPdew_psi(sel);  up = dp.ConfUpper_psi(sel) - dp.BulkPdew_psi(sel);
    th = R(R.Fluid_ID == fid(i) & R.RadiusSet == "thesis", :);
    plot(s.r_nm, s.Shift_model_psi, '-', 'Color', [0.05 0.49 0.48], 'LineWidth', 1.6, ...
         'DisplayName', "model " + cfg.preset);
    plot(th.r_nm, th.Shift_model_psi, 'o', 'Color', [0.05 0.49 0.48], 'MarkerFaceColor', 'w', ...
         'MarkerSize', 6, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    errorbar(dp.PoreRadius_nm(sel), 0.5 * (lo + up), 0.5 * (up - lo), 'ks', 'MarkerFaceColor', 'k', ...
             'LineWidth', 1.2, 'DisplayName', 'measured (r = 54 / 45 / 38 nm)');
    for n = 1:numel(sel)
        text(dp.PoreRadius_nm(sel(n)), lo(n), sprintf('%d psi', ns(sel(n))), 'FontSize', 8, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', [0.3 0.3 0.3]);
    end
    xlim(cfg.plotR_nm);
    yl = [min([lo; s.Shift_model_psi]), max([up; s.Shift_model_psi])];
    ylim([max(0, yl(1) - 0.15 * diff(yl)), yl(2) + 0.05 * diff(yl)]);
    xlabel('Mean pore radius [nm]');  ylabel('P_{dew,conf} - P_{dew,bulk} [psi]');
    title(label(i) + " (" + fid(i) + "), EF2");
    legend('Location', 'northeast', 'FontSize', 8, 'Box', 'off');
end
exportgraphics(gcf, fullfile(cfg.outDir, "netstress_test.png"), 'Resolution', 200);
fprintf('Written: %s\n', cfg.outDir);

%% Local functions -------------------------------------------------------------------------
function r_nm = localBackCalc(runner, c, K, Ptarget, bracket_nm, key)
    % Pore radius [nm] at which the model confined Pdew equals Ptarget [Pa]. The model
    % Pdew falls monotonically with r. A log-spaced scan from large to small r (warm
    % start follows the rising shift; failed solves at very small r are skipped) finds
    % the first sign change, then fzero refines it in log(r). NaN when no sign change
    % is found inside the bracket (e.g. target below the model bulk Pdew).
    g = @(lr) localPdew(runner, c, K, exp(lr), key) - Ptarget;
    grid = log(logspace(log10(bracket_nm(2)), log10(bracket_nm(1)), 16));
    r_nm = NaN;  prev = [];
    for lr = grid
        v = g(lr);
        if ~isfinite(v), continue; end
        if ~isempty(prev) && sign(v) ~= sign(prev(2))
            lrr = fzero(g, sort([prev(1), lr]), optimset('TolX', 1e-4));
            r_nm = exp(lrr);
            return
        end
        prev = [lr, v];
    end
end

function P = localPdew(runner, c, K, r_nm, key)
    res = runner.solve(c, "FWIOnly", K, 'r', r_nm * 1e-9, 'Key', key);
    P = NaN;  if res.ok, P = res.P; end
end
