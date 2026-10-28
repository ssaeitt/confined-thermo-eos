%% MLDO_CHECK  Bulk dew point AND maximum liquid dropout (MLDO) of the thesis fluids
%
% Purpose: test the bulk EOS against a second, independent measurement (MLDO,
% Doryani thesis Table 3.4, and F7 at 40/60 degC) so that the composition-adjustment
% vs BIP question can be judged on more than the bulk Pdew alone.
%
% For every row of sheet Bulk_PVT, four bulk (r_p = Inf) variants:
%   composition : "nominal"  (as prepared, Bulk_PVT)  |  "adjusted" (Mix_Props, if it differs)
%   BIP set     : "workbook" (sheet BIP)              |  "literature" (sheet BIP_Literature)
% For each: bulk Pdew (SaturationSolver, upper dew branch checked by flash), then a
% CCE from Pdew down to cfg.PminFrac * Pdew, MLDO in vol% of V_sat and its pressure.
%
% Section 5 (cfg.bipSweep): one BIP pair is varied at the NOMINAL compositions of
% every 20 degC fluid that contains the pair (default C1-nC5: F1-F4 and F10), all
% other BIPs from sheet BIP. Answers whether the composition adjustment or the BIP
% is the right handle for the bulk Pdew.
%
% Bulk only: no confinement term is involved (c_i = 0 at r_p = Inf), so k, p,
% lambda and the TIPs play no role here. Nothing is fitted.
%
% Needs (local): config/MixtureData - Complete Mixtures.xlsx with sheets
% Bulk_PVT, Mix_Props (Fluid_ID column) and BIP_Literature.

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
% --- Model variant -----------------------------------------------------------------
%   "frozen" : current calibration (workbook "MixtureData - Complete Mixtures.xlsx", preset jL)
%   "recal"  : C1-nC5 recalibration (workbook "MixtureData - Recalibrated.xlsx", preset jL2)
%   "psia"   : absolute pressures, BIPs and group TIPs refitted (workbook "MixtureData -
%              Recalibrated psia.xlsx", preset jL3; written by Recalibrate_psia.m)
cfg.variant = "psia";
variants = struct( ...
    'frozen', struct('file', "MixtureData - Complete Mixtures.xlsx",  'preset', "jL",  'basis', "recorded"), ...
    'recal',  struct('file', "MixtureData - Recalibrated.xlsx",      'preset', "jL2", 'basis', "recorded"), ...
    'psia',   struct('file', "MixtureData - Recalibrated psia.xlsx", 'preset', "jL3", 'basis', "absolute"));
vsel = variants.(char(cfg.variant));
cfg.pressureBasis = vsel.basis;        % "absolute": psig rows + 1 atm (entities.PressureBasis)
cfg.mixtureFile = fullfile(scriptDir, "config", vsel.file);
cfg.outDir      = fullfile(scriptDir, "results", "mldo_check_" + cfg.variant);
cfg.compositions = ["nominal", "adjusted"];
cfg.bipSets      = ["workbook", "literature"];
cfg.bipSheets    = struct('workbook', "BIP", 'literature', "BIP_Literature");
cfg.PminFrac     = 0.10;       % CCE down to 10 % of Pdew
cfg.nCoarse      = 40;         % CCE points below 0.95 Pdew
cfg.wallRock     = "EF2";      % only needed to build the EOS object (no effect at r = Inf)
cfg.bipSweep     = struct('pair', ["C1", "nC5"], ...
                          'values', [0.012022, 0.008, 0.005, 0.002, 0.0005, 0, -0.005], ...
                          'T', 293.15);           % set 'values' to [] to skip section 5

if ~isfile(cfg.mixtureFile)
    error('mldo:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.mixtureFile);
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

%% 1. Data ---------------------------------------------------------------------------
bp = entities.PressureBasis.readSheet(cfg.mixtureFile, "Bulk_PVT", cfg.pressureBasis);
bp = bp(~ismissing(string(bp.Fluid_ID)) & string(bp.Fluid_ID) ~= "", :);
mp = readtable(cfg.mixtureFile, 'Sheet', 'Mix_Props', 'VariableNamingRule', 'preserve');
if ~all(ismember(["Fluid_ID", "Campaign"], string(mp.Properties.VariableNames)))
    error('mldo:MixProps', 'Mix_Props needs the Fluid_ID and Campaign columns.');
end
rock = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallRock);

%% 2. Runs -------------------------------------------------------------------------------
rows = {};
curves = {};
for i = 1:height(bp)
    fid  = string(bp.Fluid_ID(i));
    mix  = string(bp.Mixture(i));
    T    = bp.T(i);
    zNom = localParseZ(bp.mole_frac(i));
    zAdj = localAdjusted(mp, fid);
    meas = struct('Pdew', bp.Bulk_Pdew(i), 'MLDO', bp.MLDO_pct(i), 'P_MLDO', bp.P_MLDO_psi(i));

    for compLabel = cfg.compositions
        if compLabel == "nominal"
            z = zNom;
        else
            if isempty(zAdj) || max(abs(zAdj - zNom)) < 1e-12, continue; end
            z = zAdj;
        end
        for bipLabel = cfg.bipSets
            fl  = localFluid(cfg.mixtureFile, mix, cfg.bipSheets.(char(bipLabel)));
            eos = thermo.ConfinedEOS(fl, rock);
            sat = solvers.SaturationSolver(eos);
            pt  = solvers.PTFlash(eos);
            pvt = solvers.PVTExperiments(pt);

            [res, branch] = localUpperDew(sat, pt, T, meas.Pdew * psi, z);
            tag = sprintf('%s | %.2f C | %s | %s', fid, T - 273.15, compLabel, bipLabel);
            if ~res.converged
                rows(end+1, :) = {fid, mix, T - 273.15, compLabel, bipLabel, localZStr(z), ...
                    meas.Pdew, NaN, NaN, meas.MLDO, NaN, NaN, meas.P_MLDO, NaN, "dew: " + res.reason}; %#ok<AGROW>
                fprintf('%-45s dew point failed: %s\n', tag, res.reason);
                continue;
            end
            Pd = res.P;
            cce = pvt.cce(T, z, Inf, Pd, solvers.PVTExperiments.pressureGrid(Pd, cfg.PminFrac * Pd, cfg.nCoarse));
            [mldo, pm] = solvers.PVTExperiments.mldo(cce);
            nFail = nnz(cce.reason ~= "converged" & cce.reason ~= "stable");

            status = branch;
            if nFail > 0, status = status + sprintf(" | %d CCE flash failures", nFail); end
            rows(end+1, :) = {fid, mix, T - 273.15, compLabel, bipLabel, localZStr(z), ...
                meas.Pdew, Pd / psi, 100 * (Pd / psi / meas.Pdew - 1), meas.MLDO, mldo, mldo - meas.MLDO, ...
                meas.P_MLDO, pm, status}; %#ok<AGROW>
            curves(end+1, :) = {fid, T, compLabel, bipLabel, cce.P_psi, cce.LDO_pct_Vsat, meas}; %#ok<AGROW>
            fprintf('%-45s Pdew %8.1f (exp %6.0f, %+6.2f %%) | MLDO %5.1f %% (exp %4.1f) at %6.0f psi\n', ...
                tag, Pd / psi, meas.Pdew, 100 * (Pd / psi / meas.Pdew - 1), mldo, meas.MLDO, pm);
        end
    end
end

R = cell2table(rows, 'VariableNames', {'Fluid', 'Mixture', 'T_C', 'Composition', 'BIPset', 'z', ...
    'Pdew_exp_psi', 'Pdew_model_psi', 'Pdew_err_pct', 'MLDO_exp_pct', 'MLDO_model_pct', 'MLDO_diff_pts', ...
    'P_MLDO_exp_psi', 'P_MLDO_model_psi', 'Status'});
R = convertvars(R, {'Fluid', 'Mixture', 'Composition', 'BIPset', 'z', 'Status'}, 'string');
disp(R);
writetable(R, fullfile(cfg.outDir, "mldo_check.csv"));

%% 3. Summary per variant (20 degC fluids) -------------------------------------------------
fprintf('\nSummary at 20 degC (mean absolute errors over the fluids each variant covers):\n');
S = {};
for compLabel = cfg.compositions
    for bipLabel = cfg.bipSets
        sel = R.Composition == compLabel & R.BIPset == bipLabel & abs(R.T_C - 20) < 1e-6 ...
              & isfinite(R.Pdew_model_psi);
        if ~any(sel), continue; end
        clean = sel & ~contains(R.Status, "failure");    % MLDO only from complete CCEs
        S(end+1, :) = {compLabel, bipLabel, nnz(sel), mean(abs(R.Pdew_err_pct(sel))), ...
            mean(abs(R.MLDO_diff_pts(clean)), 'omitnan'), mean(R.MLDO_diff_pts(clean), 'omitnan'), ...
            nnz(sel & ~clean)}; %#ok<AGROW>
    end
end
S = cell2table(S, 'VariableNames', {'Composition', 'BIPset', 'nFluids', 'MAE_Pdew_pct', ...
    'MAE_MLDO_pts', 'Bias_MLDO_pts', 'nExcluded_CCE_failures'});
S = convertvars(S, {'Composition', 'BIPset'}, 'string');
disp(S);
writetable(S, fullfile(cfg.outDir, "mldo_summary.csv"));

%% 4. Figures --------------------------------------------------------------------------------
fids = unique(string(curves(:, 1)), 'stable');
nF = numel(fids);
nCol = 4;  nRow = ceil(nF / nCol);
figure('Name', 'Bulk CCE liquid dropout vs measured MLDO', 'Position', [50 50 1500 320 * nRow]);
tiledlayout(nRow, nCol, 'TileSpacing', 'compact');
styles = struct('nominal_workbook', '-', 'nominal_literature', '--', ...
                'adjusted_workbook', '-.', 'adjusted_literature', ':');
for f = 1:nF
    nexttile; hold on; box on;
    sel = find(string(curves(:, 1)) == fids(f) & abs(cell2mat(curves(:, 2)) - 293.15) < 1e-6).';
    if isempty(sel), sel = find(string(curves(:, 1)) == fids(f)).'; end
    for k = sel
        key = curves{k, 3} + "_" + curves{k, 4};
        plot(curves{k, 5}, curves{k, 6}, styles.(char(key)), 'LineWidth', 1.3, ...
             'DisplayName', curves{k, 3} + ", " + curves{k, 4} + " BIP");
    end
    m = curves{sel(1), 7};
    if isfinite(m.P_MLDO)
        plot(m.P_MLDO, m.MLDO, 'kp', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'DisplayName', 'measured MLDO');
    else
        yline(m.MLDO, 'k:', 'DisplayName', 'measured MLDO');
    end
    xline(m.Pdew, 'k--', 'HandleVisibility', 'off');
    title(fids(f));
    xlabel('P [psi]'); ylabel('LDO [% of V_{sat}]');
    if f == 1, legend('Location', 'northwest', 'FontSize', 7); end
end
exportgraphics(gcf, fullfile(cfg.outDir, "mldo_check.png"), 'Resolution', 200);
fprintf('\nWritten: %s\n', cfg.outDir);

%% 5. BIP sweep at nominal compositions ----------------------------------------------------------
if ~isempty(cfg.bipSweep.values)
    pr = cfg.bipSweep.pair;
    fprintf('\nBIP sweep %s-%s at nominal compositions (all other BIPs from sheet BIP)\n', pr(1), pr(2));
    swRows = {};
    for i = 1:height(bp)
        mix = string(bp.Mixture(i));
        if abs(bp.T(i) - cfg.bipSweep.T) > 1e-6 || ~all(ismember(pr, split(mix, "-"))), continue; end
        z = localParseZ(bp.mole_frac(i));
        fl0 = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile, mix);
        for v = cfg.bipSweep.values
            fl  = localSetBIP(fl0, pr, v);
            eos = thermo.ConfinedEOS(fl, rock);
            sat = solvers.SaturationSolver(eos);
            pt  = solvers.PTFlash(eos);
            [res, branch] = localUpperDew(sat, pt, bp.T(i), bp.Bulk_Pdew(i) * psi, z);
            [mldo, pm, nFail] = deal(NaN, NaN, 0);
            if res.converged
                pvt = solvers.PVTExperiments(pt);
                cce = pvt.cce(bp.T(i), z, Inf, res.P, ...
                      solvers.PVTExperiments.pressureGrid(res.P, cfg.PminFrac * res.P, cfg.nCoarse));
                [mldo, pm] = solvers.PVTExperiments.mldo(cce);
                nFail = nnz(cce.reason ~= "converged" & cce.reason ~= "stable");
            end
            Pm = res.P / psi;  if ~res.converged, Pm = NaN; end
            swRows(end+1, :) = {string(bp.Fluid_ID(i)), mix, v, bp.Bulk_Pdew(i), Pm, ...
                100 * (Pm / bp.Bulk_Pdew(i) - 1), bp.MLDO_pct(i), mldo, mldo - bp.MLDO_pct(i), ...
                bp.P_MLDO_psi(i), pm, nFail, branch}; %#ok<AGROW>
            fprintf('  %-4s BIP %+.4f : Pdew %+6.2f %% | MLDO %5.1f (exp %4.1f)\n', ...
                    string(bp.Fluid_ID(i)), v, 100 * (Pm / bp.Bulk_Pdew(i) - 1), mldo, bp.MLDO_pct(i));
        end
    end
    W = cell2table(swRows, 'VariableNames', {'Fluid', 'Mixture', 'BIP', 'Pdew_exp_psi', 'Pdew_model_psi', ...
        'Pdew_err_pct', 'MLDO_exp_pct', 'MLDO_model_pct', 'MLDO_diff_pts', 'P_MLDO_exp_psi', ...
        'P_MLDO_model_psi', 'nCCEFailures', 'Branch'});
    W = convertvars(W, {'Fluid', 'Mixture', 'Branch'}, 'string');
    writetable(W, fullfile(cfg.outDir, "mldo_bip_sweep.csv"));

    % Per BIP value: RMS Pdew error over the binaries and over all fluids, MAE of MLDO
    vals = cfg.bipSweep.values(:);
    isBin = W.Mixture == strjoin(pr, "-");
    [rmsBin, rmsAll, maeM] = deal(nan(numel(vals), 1));
    for k = 1:numel(vals)
        s = abs(W.BIP - vals(k)) < 1e-12;
        rmsBin(k) = sqrt(mean(W.Pdew_err_pct(s & isBin).^2, 'omitnan'));
        rmsAll(k) = sqrt(mean(W.Pdew_err_pct(s).^2, 'omitnan'));
        maeM(k)   = mean(abs(W.MLDO_diff_pts(s & W.nCCEFailures == 0)), 'omitnan');
    end
    WS = table(vals, rmsBin, rmsAll, maeM, 'VariableNames', {'BIP', 'RMS_Pdew_err_binaries_pct', ...
        'RMS_Pdew_err_all_pct', 'MAE_MLDO_pts'});
    disp(WS);
    writetable(WS, fullfile(cfg.outDir, "mldo_bip_sweep_summary.csv"));

    figure('Name', 'BIP sweep', 'Position', [80 80 1100 420]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    fl_ = unique(W.Fluid, 'stable');
    nexttile; hold on; box on;
    for f = fl_.'
        s = W.Fluid == f;  plot(W.BIP(s), W.Pdew_err_pct(s), 'o-', 'DisplayName', f);
    end
    yline(0, 'k:', 'HandleVisibility', 'off');
    yline([-1 1], 'k--', 'HandleVisibility', 'off');
    xlabel("BIP " + pr(1) + "-" + pr(2)); ylabel('Bulk P_{dew} error [%]'); legend('Location', 'best');
    title('Bulk dew point (nominal compositions)');
    nexttile; hold on; box on;
    for f = fl_.'
        s = W.Fluid == f;  plot(W.BIP(s), W.MLDO_diff_pts(s), 'o-', 'DisplayName', f);
    end
    yline(0, 'k:', 'HandleVisibility', 'off');
    xlabel("BIP " + pr(1) + "-" + pr(2)); ylabel('MLDO model - measured [points]'); legend('Location', 'best');
    title('Maximum liquid dropout');
    exportgraphics(gcf, fullfile(cfg.outDir, "mldo_bip_sweep.png"), 'Resolution', 200);
end

%% Local functions ----------------------------------------------------------------------------
function fl = localSetBIP(fl, pair, value)
    K = fl.BIPMatrix;
    a = fl.componentIndex(pair(1));  b = fl.componentIndex(pair(2));
    K(a, b) = value;  K(b, a) = value;
    fl = entities.FluidProperties(fl.ComponentNames, fl.Tc, fl.Pc, fl.omega, fl.MW, fl.Vc, ...
        fl.Parachor, fl.ZRackett, fl.VolumeShift, fl.LJ_Size, fl.LJ_Energy, K);
end

function z = localParseZ(s)
    z = str2double(split(erase(string(s), ["[", "]", " "]), ","));
    z = z(:) / sum(z);
end

function s = localZStr(z)
    s = "[" + strjoin(compose("%.4g", z(:).'), ",") + "]";
end

function zAdj = localAdjusted(mp, fid)
    % Composition used in Mix_Props (early campaign) for this fluid, if any.
    sel = string(mp.Fluid_ID) == fid & string(mp.Campaign) == "early";
    zAdj = [];
    if any(sel)
        zAdj = localParseZ(mp.mole_frac(find(sel, 1)));
    end
end

function fl = localFluid(file, mixture, bipSheet)
    % FluidProperties with the BIP matrix taken from bipSheet (rows/columns by name).
    fl = entities.FluidProperties.loadFromWorkbook(file, mixture);
    if bipSheet == "BIP", return; end
    C = readcell(file, 'Sheet', bipSheet);
    colNames = string(C(1, 2:end));
    nm = fl.ComponentNames;
    K = zeros(fl.NC);
    for a = 1:fl.NC
        r = find(cellfun(@(v) (ischar(v) || isstring(v)) && strtrim(string(v)) == nm(a), C(2:end, 1)), 1) + 1;
        if isempty(r), error('mldo:BIPRow', '%s: no row for %s.', bipSheet, nm(a)); end
        for b = 1:fl.NC
            cc = find(strtrim(colNames) == nm(b), 1) + 1;
            if isempty(cc), error('mldo:BIPCol', '%s: no column for %s.', bipSheet, nm(b)); end
            K(a, b) = C{r, cc};
        end
    end
    fl = entities.FluidProperties(fl.ComponentNames, fl.Tc, fl.Pc, fl.omega, fl.MW, fl.Vc, ...
        fl.Parachor, fl.ZRackett, fl.VolumeShift, fl.LJ_Size, fl.LJ_Energy, K);
end

function [res, branch] = localUpperDew(sat, pt, T, Pguess, z)
    % Dew point on the UPPER (retrograde) branch: single phase just above, two below.
    branch = "not checked";
    for g = [1, 1.3, 0.8, 1.6]
        res = sat.solve("dew", T, g * Pguess, z, Inf);
        if ~res.converged, continue; end
        above = pt.solve(res.P * 1.002, T, z, Inf);
        below = pt.solve(res.P * 0.998, T, z, Inf, 'K_seed', res.K);
        if above.phases == 1 && below.phases == 2
            branch = "upper dew";
            return;
        end
        branch = "lower dew branch rejected";
    end
    if res.converged && branch ~= "upper dew"
        res.converged = false;
        res.reason = "only a lower-branch dew point was found";
    end
end
