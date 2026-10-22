%% RECALIBRATE_C1NC5  Re-derive the recalibrated C1-nC5 inputs in MATLAB (variant "recal")
%
% Background: the MLDO check showed that the adjusted C1-nC5 compositions reproduce the
% bulk Pdew but not the measured MLDO (F2: 17.3 % vs 25 %), whereas the NOMINAL
% compositions with a near-zero C1-nC5 BIP reproduce both. This script re-derives, in
% the MATLAB code base, every number that changes under that recalibration:
%
%   1. BIP C1-nC5   : least squares on the bulk Pdew of F1-F4 (nominal z, sheet Bulk_PVT).
%                     Reported for transparency; the workbook keeps 0.0005.
%   2. BIP nC5-nC10 : root of  Pdew_bulk,model(F10) = 4538 psi  (C1-nC5 = 0.0005).
%   3. C1-nC5 group TIP (FWI-Only, k, p, lambda frozen at jL): minimum of the RMSE of the
%      absolute confined Pdew error over F2-F4 on EF2 -- the same objective as the jL
%      fit (Fit_parameters.m, pdTarget "absolute"). The shift-based optimum is also
%      reported.
%   4. Rock-offset correlation refit with the new group TIP (thermo.RockOffsetTIP "jL2").
%
% Each MATLAB value is compared with the value written into the recalibrated workbook /
% the jL2 preset; agreement to the shown tolerance confirms the Python recalibration.
%
% Needs (local): config/MixtureData - Recalibrated.xlsx

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

cfg.file   = fullfile(scriptDir, "config", "MixtureData - Recalibrated.xlsx");
cfg.outDir = fullfile(scriptDir, "results", "recalibration");
cfg.T      = 293.15;
cfg.r_EF2  = 54e-9;
cfg.expected = struct('bipC1nC5', 0.0005, 'bipC5C10', 0.000746, ...
                      'tip', -4.891452765367, 'tipShift', -4.961798381617);
if ~isfile(cfg.file), error('recal:NoWorkbook', 'Not found: %s', cfg.file); end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi  = entities.DewPointDataset.psi2Pa;
rock = entities.RockProperties.loadFromWorkbook(cfg.file, "EF2");
bp   = readtable(cfg.file, 'Sheet', 'Bulk_PVT', 'VariableNamingRule', 'preserve');
bp   = bp(abs(bp.T - cfg.T) < 1e-6, :);

%% 1. C1-nC5 BIP from the bulk Pdew of F1-F4 (nominal compositions) -------------------
selC5 = find(string(bp.Mixture) == "C1-nC5").';
rmsC5 = @(k) sqrt(mean(arrayfun(@(i) localBulkErr(cfg, rock, bp(i, :), "C1", "nC5", k), selC5).^2));
[kC5, fC5] = fminbnd(rmsC5, -0.01, 0.02, optimset('TolX', 1e-7));
fprintf('1. BIP C1-nC5, least squares on F1-F4 bulk Pdew: %.6f (RMS %.3f %%); workbook uses %.6f (RMS %.3f %%)\n', ...
        kC5, fC5, cfg.expected.bipC1nC5, rmsC5(cfg.expected.bipC1nC5));

%% 2. nC5-nC10 BIP from the F10 bulk Pdew ------------------------------------------------
iF10 = find(string(bp.Fluid_ID) == "F10", 1);
k510 = fzero(@(k) localBulkErr(cfg, rock, bp(iF10, :), "nC5", "nC10", k), [-0.05, 0.05], ...
             optimset('TolX', 1e-9));
fprintf('2. BIP nC5-nC10 for F10 bulk Pdew = %g psi: %.6f (workbook %.6f)\n', ...
        bp.Bulk_Pdew(iF10), k510, cfg.expected.bipC5C10);

%% 3. C1-nC5 group TIP, FWI-Only, frozen k, p, lambda ------------------------------------
jL = thermo.RockOffsetTIP.frozen("jL", cfg.file);              % k, p, lambda only
dp = entities.DewPointDataset.loadFromWorkbook(cfg.file, 'DefaultT', cfg.T);
idx = dp.findCases("C1-nC5", "EF2");
runner = solvers.DewPointRunner(cfg.file, containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9}));
runner.engine("C1-nC5", "EF2");
runner.setFWI(jL.k, jL.p, jL.lambda);
errAbs   = @(t) localGroupErr(runner, dp, idx, t, cfg.r_EF2, "absolute");
errShift = @(t) localGroupErr(runner, dp, idx, t, cfg.r_EF2, "shift");
[tAbs, fAbs]     = fminbnd(@(t) sqrt(mean(errAbs(t).^2)),   -12, -1, optimset('TolX', 1e-7));
[tShift, fShift] = fminbnd(@(t) sqrt(mean(errShift(t).^2)), -12, -1, optimset('TolX', 1e-7));
fprintf('3. C1-nC5 TIP (absolute objective): %.6f, RMSE %.2f psi  (jL2 preset %.6f, diff %.1e)\n', ...
        tAbs, fAbs, cfg.expected.tip, tAbs - cfg.expected.tip);
fprintf('   C1-nC5 TIP (shift objective)   : %.6f, RMSE %.2f psi  (Python %.6f)\n', ...
        tShift, fShift, cfg.expected.tipShift);
eA = errAbs(tAbs);
for n = 1:numel(idx)
    c = dp.getCase(idx(n));
    fprintf('   %s z = [%s]: confined error %+.1f psi (%+.2f %% of the measured mid-point)\n', ...
            c.Mixture, num2str(c.z.', '%.4f '), eA(n), 100 * eA(n) / (c.PdewConfMid / psi));
end

%% 4. Rock-offset correlation (jL vs jL2) ------------------------------------------------------
t1 = thermo.RockOffsetTIP.frozen("jL",  cfg.file);
t2 = thermo.RockOffsetTIP.frozen("jL2", cfg.file);
fprintf('\n4. %s\n   %s\n', t1.summary(), t2.summary());

T = table(["BIP C1-nC5 (LSQ, F1-F4)"; "BIP nC5-nC10 (F10 bulk)"; "TIP C1-nC5 (absolute)"; ...
           "TIP C1-nC5 (shift)"; "alpha_EF2 jL2"; "alpha_B1 jL2"; "B jL2"], ...
          [kC5; k510; tAbs; tShift; t2.alphaEF2; t2.alphaB1; t2.B], ...
          [cfg.expected.bipC1nC5; cfg.expected.bipC5C10; cfg.expected.tip; cfg.expected.tipShift; ...
           NaN; NaN; NaN], ...
          'VariableNames', {'Quantity', 'MATLAB', 'Workbook_or_preset'});
disp(T);
writetable(T, fullfile(cfg.outDir, "recalibration.csv"));

%% Local functions -------------------------------------------------------------------------
function e = localBulkErr(cfg, rock, row, a, b, k)
    % Bulk Pdew error [%] of one Bulk_PVT row with BIP(a,b) = k.
    fl  = entities.FluidProperties.loadFromWorkbook(cfg.file, string(row.Mixture));
    K   = fl.BIPMatrix;
    i = fl.componentIndex(a);  j = fl.componentIndex(b);
    K(i, j) = k;  K(j, i) = k;
    fl  = entities.FluidProperties(fl.ComponentNames, fl.Tc, fl.Pc, fl.omega, fl.MW, fl.Vc, ...
          fl.Parachor, fl.ZRackett, fl.VolumeShift, fl.LJ_Size, fl.LJ_Energy, K);
    z   = str2double(split(erase(string(row.mole_frac), ["[", "]", " "]), ","));
    z   = z(:) / sum(z);
    sat = solvers.SaturationSolver(thermo.ConfinedEOS(fl, rock));
    res = sat.solve("dew", row.T, row.Bulk_Pdew * 6894.757, z, Inf);
    if ~res.converged
        e = NaN;  return;
    end
    e = 100 * (res.P / 6894.757 / row.Bulk_Pdew - 1);
end

function e = localGroupErr(runner, dp, idx, t, r, target)
    % Confined Pdew error [psi] of every case in idx with the binary TIP t.
    psi = entities.DewPointDataset.psi2Pa;
    e = nan(numel(idx), 1);
    for n = 1:numel(idx)
        c = dp.getCase(idx(n));
        res = runner.solve(c, "FWIOnly", [0, t; t, 0], 'r', r);
        if ~res.ok, e(n) = 1e3; continue; end            % failure penalty (as in Fit_parameters)
        if target == "absolute"
            e(n) = (res.P - c.PdewConfMid) / psi;
        else
            e(n) = res.shift_psi - (c.PdewConfMid - c.PdewBulk) / psi;
        end
    end
end
