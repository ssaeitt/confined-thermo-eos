%% TEMPERATURE_TEST  Blind test of the frozen confined model away from 20 degC
%
% Data: F7 (C1-nC8 [0.923, 0.077]) on EF2, net stress 500 psi, at 20/40/60 degC
% (Doryani thesis Table 4.9; sheet Mix_Props_T). Only the 20 degC point belongs to
% the calibration set; 40 and 60 degC are blind.
%
% Model: frozen FWI-Only preset (thermo.RockOffsetTIP). Nothing is refitted. The
% model's temperature dependence comes from the PR alpha(T) and the wall term
% W_i = eps (1 - exp(-eps/RT)); the TIP descriptor is T-independent for jL
% (inv_wbar) and re-evaluated at T for w1 (Tstar).
%
% Reports, per temperature:
%   bulk Pdew (model vs measured), confined Pdew (model vs measured range),
%   shift on both bases (vs model bulk, vs measured bulk), bulk/confined MLDO
%   (model) vs measured bulk MLDO.
% Plus a dense temperature sweep of the model shift, to show its trend between and
% beyond the measured temperatures (the measured shift rises then falls).

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
% --- Model variant -----------------------------------------------------------------
%   "frozen" : current calibration (workbook "MixtureData - Complete Mixtures.xlsx", preset jL)
%   "recal"  : C1-nC5 recalibration (workbook "MixtureData - Recalibrated.xlsx", preset jL2)
cfg.variant = "frozen";
variants = struct( ...
    'frozen', struct('file', "MixtureData - Complete Mixtures.xlsx", 'preset', "jL"), ...
    'recal',  struct('file', "MixtureData - Recalibrated.xlsx",     'preset', "jL2"));
vsel = variants.(char(cfg.variant));
cfg.mixtureFile = fullfile(scriptDir, "config", vsel.file);
cfg.sheet       = "Mix_Props_T";
cfg.outDir      = fullfile(scriptDir, "results", "temperature_test_" + cfg.variant);
cfg.presets     = vsel.preset;           % w1 was rejected by this test; add "w1" to repeat it
if cfg.variant == "frozen", cfg.presets = ["jL", "w1"]; end
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
cfg.sweep_C     = 10:5:80;               % model-only temperature sweep [degC]
cfg.mldo        = true;                  % CCE at each measured temperature
cfg.PminFrac    = 0.10;

if ~isfile(cfg.mixtureFile)
    error('ttest:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.mixtureFile);
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

dpT = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'Sheet', cfg.sheet);
mt  = readtable(cfg.mixtureFile, 'Sheet', cfg.sheet, 'VariableNamingRule', 'preserve');
fprintf('%d rows in %s at T = %s K\n\n', dpT.NumCases, cfg.sheet, mat2str(dpT.T.'));

%% 1. Measured temperatures ------------------------------------------------------------
rows = {};
sweepRows = {};
for preset = cfg.presets
    tip = thermo.RockOffsetTIP.frozen(preset, cfg.mixtureFile);
    fprintf('%s\n', tip.summary());
    runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius, 'maxShift_psi', Inf);

    for i = 1:dpT.NumCases
        c   = dpT.getCase(i);
        r   = cfg.poreRadius(char(c.Rock));
        eng = runner.engine(c.Mixture, c.Rock);
        runner.setFWI(tip.k, tip.p, tip.lambda);
        K   = tip.matrix(eng.eos.Fluid, c.Rock, c.T);
        res = runner.solve(c, "FWIOnly", K, 'r', r);

        exp_ = struct('bulk', dpT.BulkPdew_psi(i), 'lo', dpT.ConfLower_psi(i), 'up', dpT.ConfUpper_psi(i));
        [mB, mC, pB, pC] = deal(NaN);
        if res.ok && cfg.mldo
            pvt = solvers.PVTExperiments(solvers.PTFlash(eng.eos));
            Pb  = res.PdewBulkModel;
            cB  = pvt.cce(c.T, c.z, Inf, Pb, solvers.PVTExperiments.pressureGrid(Pb, cfg.PminFrac * Pb, 40));
            cC  = pvt.cce(c.T, c.z, r, res.P, solvers.PVTExperiments.pressureGrid(res.P, cfg.PminFrac * res.P, 40));
            [mB, pB] = solvers.PVTExperiments.mldo(cB);
            [mC, pC] = solvers.PVTExperiments.mldo(cC);
        end
        Pc = res.P / psi;  Pbm = res.PdewBulkModel / psi;
        inRange = Pc >= exp_.lo && Pc <= exp_.up;
        shM = Pc - Pbm;                                  % vs model bulk
        inShift = shM >= exp_.lo - exp_.bulk && shM <= exp_.up - exp_.bulk;
        role = "blind";  if abs(c.T - 293.15) < 1e-6, role = "calibration"; end
        rows(end+1, :) = {preset, c.T - 273.15, role, exp_.bulk, Pbm, 100 * (Pbm / exp_.bulk - 1), ...
            exp_.lo, exp_.up, Pc, 100 * (Pc / (0.5 * (exp_.lo + exp_.up)) - 1), inRange, ...
            exp_.lo - exp_.bulk, exp_.up - exp_.bulk, shM, Pc - exp_.bulk, inShift, ...
            mt.MLDO_pct(i), mB, pB, mC, pC, string(res.reason)}; %#ok<AGROW>
        fprintf(['  %4.0f C (%s): bulk %7.1f (exp %4.0f) | conf %7.1f (exp %4.0f-%4.0f) | shift %+6.1f ' ...
                 '(exp %+4.0f to %+4.0f) | MLDO bulk %4.1f / conf %4.1f (exp bulk %2.0f)\n'], ...
                c.T - 273.15, role, Pbm, exp_.bulk, Pc, exp_.lo, exp_.up, shM, ...
                exp_.lo - exp_.bulk, exp_.up - exp_.bulk, mB, mC, mt.MLDO_pct(i));
    end

    % Model-only sweep: same fluid and rock, warm-started from the previous T
    c0 = dpT.getCase(1);
    for Tc = cfg.sweep_C
        c = c0;  c.T = Tc + 273.15;
        eng = runner.engine(c.Mixture, c.Rock);
        K   = tip.matrix(eng.eos.Fluid, c.Rock, c.T);
        res = runner.solve(c, "FWIOnly", K, 'r', cfg.poreRadius(char(c.Rock)), 'Key', "sweep|" + preset);
        sweepRows(end+1, :) = {preset, Tc, res.PdewBulkModel / psi, res.P / psi, ...
                               (res.P - res.PdewBulkModel) / psi, string(res.reason)}; %#ok<AGROW>
    end
    fprintf('\n');
end

R = cell2table(rows, 'VariableNames', {'Preset', 'T_C', 'Role', 'Bulk_exp_psi', 'Bulk_model_psi', ...
    'Bulk_err_pct', 'Conf_exp_lo_psi', 'Conf_exp_up_psi', 'Conf_model_psi', 'Conf_err_vs_mid_pct', ...
    'Conf_in_range', 'Shift_exp_lo_psi', 'Shift_exp_up_psi', 'Shift_model_vs_modelBulk_psi', ...
    'Shift_model_vs_expBulk_psi', 'Shift_in_range', 'MLDO_exp_bulk_pct', 'MLDO_model_bulk_pct', ...
    'P_MLDO_model_bulk_psi', 'MLDO_model_conf_pct', 'P_MLDO_model_conf_psi', 'Reason'});
R = convertvars(R, {'Preset', 'Role', 'Reason'}, 'string');
Sw = cell2table(sweepRows, 'VariableNames', {'Preset', 'T_C', 'Bulk_model_psi', 'Conf_model_psi', ...
    'Shift_model_psi', 'Reason'});
Sw = convertvars(Sw, {'Preset', 'Reason'}, 'string');
disp(R);
writetable(R,  fullfile(cfg.outDir, "temperature_test.csv"));
writetable(Sw, fullfile(cfg.outDir, "temperature_sweep.csv"));

%% 2. Figure: shift vs temperature -------------------------------------------------------------
figure('Name', 'Confined Pdew shift vs temperature', 'Position', [80 80 1100 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on; box on;
Te = dpT.T - 273.15;
lo = dpT.ConfLower_psi - dpT.BulkPdew_psi;  up = dpT.ConfUpper_psi - dpT.BulkPdew_psi;
errorbar(Te, 0.5 * (lo + up), 0.5 * (up - lo), 'ks', 'MarkerFaceColor', 'k', 'LineWidth', 1.2, ...
         'DisplayName', 'measured (EF2, F7)');
for preset = cfg.presets
    s = Sw(Sw.Preset == preset, :);
    plot(s.T_C, s.Shift_model_psi, '-', 'LineWidth', 1.4, 'DisplayName', "model " + preset);
end
xline(20, ':', 'calibration T', 'HandleVisibility', 'off');
xlabel('T [\circC]'); ylabel('P_{dew,conf} - P_{dew,bulk} [psi]'); legend('Location', 'best');
title('Shift (model vs model bulk; measured vs measured bulk)');
nexttile; hold on; box on;
plot(Te, dpT.BulkPdew_psi, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'measured bulk');
errorbar(Te, 0.5 * (dpT.ConfLower_psi + dpT.ConfUpper_psi), 0.5 * (dpT.ConfUpper_psi - dpT.ConfLower_psi), ...
         'ks', 'DisplayName', 'measured confined');
for preset = cfg.presets
    s = Sw(Sw.Preset == preset, :);
    plot(s.T_C, s.Bulk_model_psi, '--', 'DisplayName', "model bulk " + preset);
    plot(s.T_C, s.Conf_model_psi, '-',  'DisplayName', "model confined " + preset);
end
xlabel('T [\circC]'); ylabel('P_{dew} [psi]'); legend('Location', 'best', 'FontSize', 7);
title('Absolute dew points');
exportgraphics(gcf, fullfile(cfg.outDir, "temperature_test.png"), 'Resolution', 200);
fprintf('Written: %s\n', cfg.outDir);
