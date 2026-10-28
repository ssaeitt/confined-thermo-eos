%% RECALIBRATE_PSIA  Absolute-pressure recalibration (variant "psia", preset jL3)
%
% Laboratory answer (28 Oct 2026): the early campaign (Doryani thesis; LabView VI or
% Druck display unit) recorded gauge pressures (psig); the late campaign (Quartzdyne
% with Q-Com, from ~2017) records absolute pressures (psia). On the absolute basis the
% two C1-nC8 [0.9, 0.1] bulk values are 4124.7 and 4126 psia (1.3 psi apart).
%
% This script
%   0. copies "MixtureData - Recalibrated.xlsx" to "MixtureData - Recalibrated psia.xlsx"
%      and writes the P_basis flag column (tools/add_pressure_basis.m). Recorded values
%      are not changed; entities.PressureBasis adds 1 atm to the psig rows on loading.
%   1. refits the bulk BIPs (sheet Bulk_PVT, 20 degC, nominal compositions), each by
%      least squares on the bulk Pdew [psi], in this order:
%        C1-nC5 (F1-F4), C1-nC8 (F7-F9), C1-nC10 (F5-F6), then nC5-nC10 from F10.
%      Done on both bases ("recorded" and "absolute") so the effect of the conversion
%      alone is visible. The late B1 C8(10) bulk value is reported, not fitted.
%   2. refits the four group TIPs (FWI-Only; k, p, lambda of jL), same objective as
%      the jL / jL2 fits: RMSE of the absolute confined Pdew over the early-campaign
%      rows of the group (Recalibrate_C1nC5.m, step 3).
%        run "check" : recorded pressures, workbook BIPs  -> should give the jL2 TIPs
%        run "psia"  : absolute pressures, BIPs of step 1 -> preset jL3
%   3. fits the rock-offset correlation (LL-inv_wbar) to the jL3 group TIPs.
%   4. back-calculates the late B1 pore radius: model confined Pdew of C1-nC8
%      [0.9, 0.1] on B1 = mid-point of the measured range (absolute target, as
%      Baclcalc_rp_B1_late.m), with the jL3 correlation and BIPs.
%   5. writes into the psia workbook: sheet BIP (step-1 absolute BIPs), TIP_Groups
%      (read by thermo.RockOffsetTIP.frozen("jL3", ...)) and Calibration_Info (late
%      radius, read by Uncertainty_BIP_TIP.m), and a summary to results/recalibration_psia.
%
% Afterwards: Uncertainty_BIP_TIP.m, PVT_verification.m, MLDO_check.m,
% Temperature_test.m and NetStress_test.m with cfg.variant = "psia".
%
% Needs (local): config/MixtureData - Recalibrated.xlsx

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));
addpath(fullfile(scriptDir, "tools"));

%% 0. Configuration and workbook ---------------------------------------------------------
cfg.source = fullfile(scriptDir, "config", "MixtureData - Recalibrated.xlsx");
cfg.target = fullfile(scriptDir, "config", "MixtureData - Recalibrated psia.xlsx");
cfg.outDir = fullfile(scriptDir, "results", "recalibration_psia");
cfg.T      = 293.15;
cfg.poreRadius = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
% BIP fits, in order (each uses the BIPs fitted before it)
cfg.bipFits = struct('a', {"C1", "C1", "C1"}, 'b', {"nC5", "nC8", "nC10"}, ...
                     'fluids', {["F1", "F2", "F3", "F4"], ["F7", "F8", "F9"], ["F5", "F6"]}, ...
                     'range', {[-0.03, 0.04], [-0.01, 0.05], [0.0, 0.08]});
cfg.ternary = struct('a', "nC5", 'b', "nC10", 'fluid', "F10", 'range', [-0.05, 0.05]);
cfg.tipSearch = struct('grid', linspace(-2, 2.6, 24), 'tolX', 1e-7);   % ln(-TIP)
cfg.late = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1], 'bracket_nm', [4, 14]);

if ~isfile(cfg.source), error('recalpsia:NoWorkbook', 'Not found: %s', cfg.source); end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;
ref = thermo.RockOffsetTIP.frozen("jL2", cfg.source);          % k, p, lambda and jL2 TIPs
kpl = [ref.k, ref.p, ref.lambda];

copyfile(cfg.source, cfg.target, 'f');
rep = add_pressure_basis(cfg.target, 'Overwrite', true);
fprintf('0. %s\n', cfg.target);
disp(rep);
rock = entities.RockProperties.loadFromWorkbook(cfg.target, "EF2");    % no effect at r = Inf
FL = containers.Map();                     % fluids with the source-workbook BIPs, by mixture
for mix = ["C1-nC5", "C1-nC8", "C1-nC10", "C1-nC5-nC10"]
    FL(char(mix)) = entities.FluidProperties.loadFromWorkbook(cfg.source, mix);
end

%% 1. Bulk BIPs on both bases --------------------------------------------------------------
bases = ["recorded", "absolute"];
BIP = struct();  bipRows = {};  residRows = {};
for b = bases
    bp = entities.PressureBasis.readSheet(cfg.target, "Bulk_PVT", b);
    bp = bp(abs(bp.T - cfg.T) < 1e-6, :);
    K = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for f = cfg.bipFits
        sel = find(ismember(string(bp.Fluid_ID), f.fluids)).';
        rmsFun = @(k) sqrt(mean(arrayfun(@(i) localBulkErr(FL, rock, bp(i, :), localWith(K, f.a, f.b, k)), sel).^2));
        [k, fv] = fminbnd(rmsFun, f.range(1), f.range(2), optimset('TolX', 1e-7));
        K(char(f.a + "|" + f.b)) = k;
        wb = localWorkbookBIP(cfg.source, f.a, f.b);
        bipRows(end+1, :) = {b, f.a + "-" + f.b, strjoin(f.fluids, ","), k, fv, wb}; %#ok<SAGROW>
        fprintf('1. %-8s BIP %s-%s = %.6f (RMS %.2f psi; workbook %.6f)\n', b, f.a, f.b, k, fv, wb);
    end
    iT = find(string(bp.Fluid_ID) == cfg.ternary.fluid, 1);
    t = cfg.ternary;
    k = fzero(@(k) localBulkErr(FL, rock, bp(iT, :), localWith(K, t.a, t.b, k)), t.range, optimset('TolX', 1e-9));
    K(char(t.a + "|" + t.b)) = k;
    bipRows(end+1, :) = {b, t.a + "-" + t.b, t.fluid, k, 0, localWorkbookBIP(cfg.source, t.a, t.b)}; %#ok<SAGROW>
    fprintf('1. %-8s BIP %s-%s = %.6f (exact on %s)\n', b, t.a, t.b, k, t.fluid);
    for i = 1:height(bp)
        e = localBulkErr(FL, rock, bp(i, :), K);
        residRows(end+1, :) = {b, string(bp.Fluid_ID(i)), string(bp.Mixture(i)), bp.Bulk_Pdew(i), ...
                               bp.Bulk_Pdew(i) + e, e}; %#ok<SAGROW>
    end
    % late B1 C8(10) bulk value (Mix_Props), not fitted
    dpb = entities.DewPointDataset.loadFromWorkbook(cfg.target, 'DefaultT', cfg.T, 'PressureBasis', b);
    iL = localLateIndex(dpb, cfg.late);
    cL = dpb.getCase(iL);
    row = table(string(cL.Mixture), "[0.9,0.1]", cfg.T, cL.PdewBulk / psi, 'VariableNames', ...
                {'Mixture', 'mole_frac', 'T', 'Bulk_Pdew'});
    e = localBulkErr(FL, rock, row, K);
    residRows(end+1, :) = {b, "F8 late (B1)", cL.Mixture, cL.PdewBulk / psi, cL.PdewBulk / psi + e, e}; %#ok<SAGROW>
    BIP.(char(b)) = K;
end
BipTable = cell2table(bipRows, 'VariableNames', {'Basis', 'Pair', 'Fluids', 'BIP', 'RMS_psi', 'Workbook_BIP'});
BipTable = convertvars(BipTable, {'Basis', 'Pair', 'Fluids'}, 'string');
BulkResid = cell2table(residRows, 'VariableNames', {'Basis', 'Fluid', 'Mixture', 'Pdew_meas_psi', ...
                                                    'Pdew_model_psi', 'Resid_psi'});
BulkResid = convertvars(BulkResid, {'Basis', 'Fluid', 'Mixture'}, 'string');
disp(BipTable);  disp(BulkResid);

%% 2. Group TIPs: regression check (recorded, workbook BIPs) and psia (absolute, new BIPs) ----
runs = struct('name', {"check", "psia"}, 'basis', {"recorded", "absolute"}, 'useFitBIP', {false, true});
groups = thermo.RockOffsetTIP.Groups;
TIP = struct();  tipRows = {};
for rn = runs
    dp = entities.DewPointDataset.loadFromWorkbook(cfg.target, 'DefaultT', cfg.T, 'PressureBasis', rn.basis);
    runner = solvers.DewPointRunner(cfg.target, cfg.poreRadius);
    runner.setFWI(kpl(1), kpl(2), kpl(3));
    if rn.useFitBIP, localApplyBIP(runner, BIP.absolute); end
    t = nan(1, numel(groups));
    for g = 1:numel(groups)
        gp = split(groups(g), "|");
        idx = dp.findCases(gp(1), gp(2));
        idx = idx(dp.HasConfined(idx) & abs(dp.T(idx) - cfg.T) < 1e-6 & dp.Campaign(idx) ~= "late");
        f = localFitGroupTIP(runner, dp, idx, cfg);
        t(g) = f.TIP;
        tipRows(end+1, :) = {rn.name, groups(g), numel(idx), f.TIP, ref.GroupTIPs(g), f.rmse_psi, ...
                             strjoin(compose("%+.1f", f.resid_psi.'), " / "), f.atEdge}; %#ok<SAGROW>
        fprintf('2. %-5s %-11s TIP %.6f (jL2 %.6f) RMSE %.2f psi | residuals %s psi%s\n', rn.name, groups(g), ...
                f.TIP, ref.GroupTIPs(g), f.rmse_psi, strjoin(compose("%+.1f", f.resid_psi.'), " / "), ...
                localIf(f.atEdge, "  [at search edge]", ""));
    end
    TIP.(char(rn.name)) = t;
end
TipTable = cell2table(tipRows, 'VariableNames', {'Run', 'Group', 'nCases', 'TIP', 'TIP_jL2', 'RMSE_psi', ...
                                                 'Residuals_psi', 'AtSearchEdge'});
TipTable = convertvars(TipTable, {'Run', 'Group', 'Residuals_psi'}, 'string');

%% 3. Rock-offset correlation ---------------------------------------------------------------
tipChk = thermo.RockOffsetTIP.custom("check", TIP.check, "inv_wbar", cfg.target, kpl);
tip3   = thermo.RockOffsetTIP.custom("psia",  TIP.psia,  "inv_wbar", cfg.target, kpl);
fprintf('\n3. jL2   : %s\n   check : %s\n   psia  : %s\n', ref.summary(), tipChk.summary(), tip3.summary());

%% 4. Late B1 pore radius (absolute basis, jL3) ----------------------------------------------
dpA = entities.DewPointDataset.loadFromWorkbook(cfg.target, 'DefaultT', cfg.T, 'PressureBasis', "absolute");
cL = dpA.getCase(localLateIndex(dpA, cfg.late));
rr = solvers.DewPointRunner(cfg.target, cfg.poreRadius, 'maxShift_psi', Inf);
rr.setFWI(kpl(1), kpl(2), kpl(3));
localApplyBIP(rr, BIP.absolute);
flL = entities.FluidProperties.loadFromWorkbook(cfg.target, cL.Mixture);
flL = localFluidWith(flL, BIP.absolute);
KL = tip3.matrix(flL, "B1", cfg.T);
Pc = @(r_nm) localConfP(rr, cL, KL, r_nm * 1e-9);
scan = cfg.late.bracket_nm(1):0.5:cfg.late.bracket_nm(2);
v = arrayfun(@(r) Pc(r) - cL.PdewConfMid, scan);
rLate = NaN;  rLo = NaN;  rHi = NaN;
j = find(isfinite(v(1:end-1)) & isfinite(v(2:end)) & sign(v(1:end-1)) ~= sign(v(2:end)), 1);
if ~isempty(j)
    rLate = fzero(@(r) Pc(r) - cL.PdewConfMid, scan([j, j + 1]), optimset('TolX', 1e-4));
    lo = cL.PdewConfMid - cL.PdewConfHalfWidth;  hi = cL.PdewConfMid + cL.PdewConfHalfWidth;
    try, rHi = fzero(@(r) Pc(r) - lo, [rLate, cfg.late.bracket_nm(2)], optimset('TolX', 1e-4)); catch, end
    try, rLo = fzero(@(r) Pc(r) - hi, [cfg.late.bracket_nm(1), rLate], optimset('TolX', 1e-4)); catch, end
end
fprintf(['\n4. Late B1 C1-nC8 [0.9, 0.1]: measured %.1f-%.1f psia, bulk %.1f psia -> r_p = %.3f nm ' ...
         '(range %.3f-%.3f nm; recal value 7.827 nm)\n'], (cL.PdewConfMid - cL.PdewConfHalfWidth) / psi, ...
         (cL.PdewConfMid + cL.PdewConfHalfWidth) / psi, cL.PdewBulk / psi, rLate, rLo, rHi);
if ~isfinite(rLate), warning('recalpsia:NoRadius', 'No late-core radius in [%g, %g] nm.', cfg.late.bracket_nm); end

%% 5. Write the psia workbook and the summary -------------------------------------------------
C = readcell(cfg.target, 'Sheet', 'BIP');
miss = cellfun(@(e) isa(e, 'missing'), C);  C(miss) = {''};
hdr = string(C(1, 2:end));  lab = string(C(2:end, 1));
ks = keys(BIP.absolute);
for n = 1:numel(ks)
    pr = split(string(ks{n}), "|");
    ra = find(lab == pr(1), 1);  rb = find(lab == pr(2), 1);  ca = find(hdr == pr(1), 1);  cb = find(hdr == pr(2), 1);
    C{1 + ra, 1 + cb} = BIP.absolute(ks{n});  C{1 + rb, 1 + ca} = BIP.absolute(ks{n});
end
writecell(C, cfg.target, 'Sheet', 'BIP', 'Range', 'A1');

TG = table(groups(:), TIP.psia(:), TipTable.RMSE_psi(TipTable.Run == "psia"), ...
           TipTable.nCases(TipTable.Run == "psia"), repmat("absolute confined Pdew, FWI-Only, k/p/lambda of jL", numel(groups), 1), ...
           'VariableNames', {'Group', 'TIP', 'RMSE_psi', 'nCases', 'Objective'});
writetable(TG, cfg.target, 'Sheet', 'TIP_Groups');
info = { ...
    "PressureBasis",    "absolute (psig rows + " + sprintf('%.6f', entities.PressureBasis.Atm_psi) + " psi)"; ...
    "Source",           "Recalibrate_psia.m from MixtureData - Recalibrated.xlsx, " + string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')); ...
    "k_p_lambda",       sprintf('%.12g, %.12g, %.12g (jL)', kpl); ...
    "alpha_EF2",        sprintf('%.9f', tip3.alphaEF2); ...
    "alpha_B1",         sprintf('%.9f', tip3.alphaB1); ...
    "B",                sprintf('%.9f', tip3.B); ...
    "LateB1_rp_nm",     sprintf('%.6f', rLate); ...
    "LateB1_rp_lo_nm",  sprintf('%.6f', rLo); ...
    "LateB1_rp_hi_nm",  sprintf('%.6f', rHi)};
info = cellfun(@string, info, 'UniformOutput', false);
writetable(cell2table(info, 'VariableNames', {'Item', 'Value'}), cfg.target, 'Sheet', 'Calibration_Info');

xls = fullfile(cfg.outDir, "recalibration_psia.xlsx");
if isfile(xls), delete(xls); end
writetable(rep, xls, 'Sheet', 'PressureBasis');
writetable(BipTable, xls, 'Sheet', 'BIPs');
writetable(BulkResid, xls, 'Sheet', 'BulkResiduals');
writetable(TipTable, xls, 'Sheet', 'GroupTIPs');
Corr = table(["jL2"; "check"; "psia (jL3)"], [ref.alphaEF2; tipChk.alphaEF2; tip3.alphaEF2], ...
             [ref.alphaB1; tipChk.alphaB1; tip3.alphaB1], [ref.B; tipChk.B; tip3.B], ...
             [ref.rmseLnEF2; tipChk.rmseLnEF2; tip3.rmseLnEF2], ...
             'VariableNames', {'Model', 'alpha_EF2', 'alpha_B1', 'B', 'RMSE_ln_EF2'});
writetable(Corr, xls, 'Sheet', 'Correlation');
writetable(cell2table(info, 'VariableNames', {'Item', 'Value'}), xls, 'Sheet', 'Calibration_Info');
disp(TipTable);  disp(Corr);
fprintf('Written: %s (sheets BIP, TIP_Groups, Calibration_Info)\n         %s\n', cfg.target, xls);
chk = max(abs(TIP.check - ref.GroupTIPs) ./ abs(ref.GroupTIPs));
fprintf('Regression check (recorded pressures, workbook BIPs vs jL2 TIPs): max relative difference %.2g\n', chk);

%% Local functions ------------------------------------------------------------------------------
function K2 = localWith(K, a, b, k)
    % Copy of the BIP override map with (a, b) = k.
    K2 = containers.Map('KeyType', 'char', 'ValueType', 'any');
    ks = keys(K);
    for n = 1:numel(ks), K2(ks{n}) = K(ks{n}); end
    K2(char(a + "|" + b)) = k;
end

function k = localWorkbookBIP(file, a, b)
    fl = entities.FluidProperties.loadFromWorkbook(file, a + "-" + b);
    k = fl.BIPMatrix(fl.componentIndex(a), fl.componentIndex(b));
end

function fl = localFluidWith(fl, K)
    % FluidProperties with the BIP overrides of map K ("A|B" -> value) applied.
    B = fl.BIPMatrix;  ks = keys(K);
    for n = 1:numel(ks)
        pr = split(string(ks{n}), "|");
        if ~all(ismember(pr, string(fl.ComponentNames))), continue; end
        i = fl.componentIndex(pr(1));  j = fl.componentIndex(pr(2));
        B(i, j) = K(ks{n});  B(j, i) = B(i, j);
    end
    fl = entities.FluidProperties(fl.ComponentNames, fl.Tc, fl.Pc, fl.omega, fl.MW, fl.Vc, ...
         fl.Parachor, fl.ZRackett, fl.VolumeShift, fl.LJ_Size, fl.LJ_Energy, B);
end

function localApplyBIP(runner, K)
    ks = keys(K);
    for n = 1:numel(ks)
        pr = split(string(ks{n}), "|");
        runner.setBIP(pr(1), pr(2), K(ks{n}));
    end
end

function e = localBulkErr(FL, rock, row, K)
    % Model - measured bulk Pdew [psi] of one row (Mixture, mole_frac, T, Bulk_Pdew),
    % upper dew branch; FL: fluids by mixture (source BIPs), K: BIP overrides.
    fl  = localFluidWith(FL(char(string(row.Mixture))), K);
    z   = str2double(split(erase(string(row.mole_frac), ["[", "]", " "]), ","));
    z   = z(:) / sum(z);
    eos = thermo.ConfinedEOS(fl, rock);
    res = solvers.UncertaintyBudget.upperDew(solvers.SaturationSolver(eos), solvers.PTFlash(eos), ...
                                            row.T, row.Bulk_Pdew * 6894.757, z);
    e = NaN;
    if res.converged, e = res.P / 6894.757 - row.Bulk_Pdew; end
end

function i = localLateIndex(dp, L)
    i = find(dp.Mixture == L.Mixture & dp.Rock == L.Rock & dp.Campaign == "late" & ...
             arrayfun(@(n) numel(dp.MoleFrac{n}) == numel(L.z) && max(abs(dp.MoleFrac{n} - L.z)) < 1e-9, ...
                      (1:dp.NumCases).'), 1);
    if isempty(i), error('recalpsia:NoLate', 'Late-campaign %s on %s not found in Mix_Props.', L.Mixture, L.Rock); end
end

function P = localConfP(runner, c, K, r)
    res = runner.solve(c, "FWIOnly", K, 'r', r);
    P = NaN;  if res.ok, P = res.P; end
end

function f = localFitGroupTIP(runner, dp, idx, cfg)
    % Group TIP minimising the RMSE of the absolute confined Pdew (jL / jL2 objective).
    psi = entities.DewPointDataset.psi2Pa;
    resid = @(lt) localGroupResid(runner, dp, idx, -exp(lt), psi);
    obj = @(lt) sqrt(mean(resid(lt).^2));
    g = cfg.tipSearch.grid;  v = arrayfun(obj, g);
    [~, i] = min(v);
    lo = g(max(i - 1, 1));  hi = g(min(i + 1, numel(g)));
    lt = fminbnd(obj, lo, hi, optimset('TolX', cfg.tipSearch.tolX));
    f = struct('TIP', -exp(lt), 'resid_psi', resid(lt), 'rmse_psi', obj(lt), 'atEdge', i == 1 || i == numel(g));
end

function e = localGroupResid(runner, dp, idx, t, psi)
    e = nan(numel(idx), 1);
    for n = 1:numel(idx)
        c = dp.getCase(idx(n));
        res = runner.solve(c, "FWIOnly", [0, t; t, 0]);
        if ~res.ok, e(n) = 1e3; continue; end            % failure penalty (as in Fit_parameters)
        e(n) = (res.P - c.PdewConfMid) / psi;
    end
end

function s = localIf(tf, a, b)
    if tf, s = a; else, s = b; end
end
