%% CO2_PIPELINE  Second batch: CO2-HC mixtures with the FROZEN confined model
%
% Nothing is fitted here. The model is frozen from the HC batch:
%   primary     joint-Linearization (k, p, lambda) + geomean log-linear TIP
%   sensitivity w = 1 point of the exact-objective Pareto front
%   TIP_ij = -exp( alpha + beta / sqrt(M_i M_j) + gamma * sigma_bar_ij / r_p )
%
% Stages, in the order they must be read:
%   A  blind CO2 critical-shift test  - validates c_CO2 before any mixture
%   B  late-core r_p from the B1 C1-nC8 10% point, carried as a band
%   C  bulk CO2 dew points with LITERATURE BIPs - residual reported, not tuned
%   D  confined CO2 predictions, both frozen models x {nominal r_p, derived r_p}
%   E  CO2-C1 TIP sensitivity: raw extrapolation vs M capped at the calibration
%      edge (M_geomean(CO2-C1) = 26.6 lies below the calibrated range)
%
% The frozen parameters and coefficients are written to results/frozen_model.*
% with a tag, so every later run can be traced to exactly this freeze.

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
            error('co2:PackageNotFound', 'Folder "+%s" not found under %s.', pk, scriptDir);
        end
        addpath(fileparts(cand(1)));
    end
end

%% 0. Frozen model and configuration -----------------------------------------
cfg.configDir   = fullfile(scriptDir, "config");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.mixtureFile = fullfile(cfg.configDir, "MixtureData - Complete Mixtures.xlsx");
cfg.tcShiftFile = fullfile(cfg.configDir, "Extended_Dataset - Original.xlsx");
cfg.tcSheets    = ["All Data (All comp)", "All Data"];   % first sheet that has CO2 rows is used
cfg.outFile     = "co2_pipeline";
cfg.freezeTag   = "HC-batch-frozen-" + string(datetime('now', 'Format', 'yyyyMMdd'));

cfg.T           = 293.15;
cfg.modes       = ["FWIOnly", "Combined"];
cfg.wallRock    = "EF2";        % row used only to read mineral_E
cfg.wallMineral = "Quartz";     % calibration wall: dTc only
cfg.poreRadiusNominal = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});

% --- FROZEN models (k, p, lambda and the geomean log-linear coefficients) ---
frozen = struct( ...
    'name',   {"joint-Linearization",           "pareto-exact-w1"}, ...
    'role',   {"primary",                       "sensitivity"}, ...
    'k',      {112.42634826677,                 166.410385934507}, ...
    'p',      {1.3558179290678,                 1.5488421885245}, ...
    'lambda', {2.0916235090532,                 4.47908205631447}, ...
    'alpha',  {-2.3398972497412,                -1.89471289979616}, ...
    'beta',   {141.538993399745,                137.987879389964}, ...
    'gamma',  {-25.7758596544873,               -24.2342345026821});

% --- Late-core r_p (stage B) -------------------------------------------------
cfg.rp.case      = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.rp.bracket_nm = [4, 20];        % search range for the derived r_p
cfg.rp.band_nm    = [6, 11.25];     % band carried through the predictions
cfg.rp.mode       = "FWIOnly";

% --- Composition policy (stage C/D). STATE THE PROJECT POLICY HERE ----------
cfg.composition.normalise = true;
cfg.composition.policy = "Reported mole fractions are used as given and renormalised to unity; " + ...
    "no composition is adjusted to improve the bulk dew point.";

% --- CO2-C1 TIP sensitivity (stage E) ---------------------------------------
cfg.tip.capAtCalibrationEdge = true;   % report BOTH raw and capped

cfg.accept.maxPdewErr_pct = 1.0;

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
psi = entities.DewPointDataset.psi2Pa;

fprintf('Frozen model tag: %s\n', cfg.freezeTag);
for f = frozen
    fprintf('  %-20s (%s): k = %.6f, p = %.6f, lambda = %.6f | alpha %+.6f, beta %+.6f, gamma %+.6f\n', ...
            f.name, f.role, f.k, f.p, f.lambda, f.alpha, f.beta, f.gamma);
end
fprintf('Composition policy: %s\n\n', cfg.composition.policy);

%% 1. Data and engines ---------------------------------------------------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
fprintf('Mix_Props rows: %d (%d with confined data). Mixtures: %s\n', dp.NumCases, ...
        nnz(dp.HasConfined), strjoin(unique(dp.Mixture, 'stable'), ", "));
isCO2 = contains(dp.Mixture, "CO2");
fprintf('CO2 rows: %d\n', nnz(isCO2));
if ~any(isCO2)
    error('co2:NoCO2Rows', 'No CO2 mixtures in %s. Add the second-batch rows first.', cfg.mixtureFile);
end

runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadiusNominal);
for i = find(dp.HasConfined).'
    runner.engine(dp.Mixture(i), dp.Rock(i));
end

% Calibration-edge M for the TIP cap: smallest M_geomean among the fitted HC groups
calPairs = ["C1-nC5", "C1-nC8", "C1-nC10"];
Mcal = zeros(numel(calPairs), 1);
for i = 1:numel(calPairs)
    fl = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile, calPairs(i));
    Mcal(i) = sqrt(prod(fl.MW));
end
cfg.tip.Mcap = min(Mcal);
fprintf('TIP calibration range in M_geomean: %.2f to %.2f (cap = %.2f)\n\n', min(Mcal), max(Mcal), cfg.tip.Mcap);

%% 2. Stage A: blind CO2 critical-shift test -----------------------------------
fluidAll = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile);
rockRef  = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallRock);
q        = rockRef.mineral(cfg.wallMineral);
calRock  = entities.RockProperties("CalWall-" + q.Name, 0, q.Name, q.Epsilon_K, q.GrainDensity, 1);
eosCal   = thermo.ConfinedEOS(fluidAll, calRock);

tcCO2 = [];
for sh = cfg.tcSheets
    try
        T = entities.TcShiftDataset.loadFromWorkbook(cfg.tcShiftFile, 'Sheet', sh, 'Components', "CO2");
    catch
        continue;
    end
    if T.NumPoints > 0, tcCO2 = T; cfg.tcSheetUsed = sh; break; end
end
tcRows = {};
if isempty(tcCO2)
    warning('co2:NoCO2Tc', 'No CO2 critical-shift points found in %s. Stage A skipped.', cfg.tcShiftFile);
else
    [idx, r] = tcCO2.modelInputs(fluidAll);
    fprintf('=== Stage A: blind CO2 critical-shift test (%d points, sheet "%s") ===\n', ...
            tcCO2.NumPoints, cfg.tcSheetUsed);
    for f = frozen
        eosCal.k = f.k;  eosCal.pT_wall = f.p;  eosCal.lambda = f.lambda;
        for n = 1:tcCO2.NumPoints
            i = idx(n);
            pred = NaN;
            try
                pred = 1 - eosCal.pureCriticalPoint(i, r(n)) / eosCal.Fluid.Tc(i);
            catch ME
                if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
            end
            tcRows(end+1, :) = {f.name, f.role, tcCO2.Component(n), tcCO2.dH(n), tcCO2.Reference(n), ...
                tcCO2.dTc(n), pred, pred - tcCO2.dTc(n), ...
                100 * (pred - tcCO2.dTc(n)) / tcCO2.dTc(n)}; %#ok<AGROW>
        end
    end
    tcTable = cell2table(tcRows, 'VariableNames', {'Model', 'Role', 'Component', 'dH', 'Reference', ...
        'dTc_exp', 'dTc_pred', 'Residual', 'Error_pct'});
    tcTable = convertvars(tcTable, {'Model', 'Role', 'Component', 'Reference'}, 'string');
    disp(tcTable);
    fprintf('Stage A is a blind test of c_CO2: no CO2 information entered the frozen parameters.\n\n');
end

%% 3. Stage B: late-core r_p from the B1 C1-nC8 10% point ----------------------
iRp = dp.findCases(cfg.rp.case.Mixture, cfg.rp.case.Rock, cfg.rp.case.z);
rpRows = {};
rpDerived = containers.Map();
if isempty(iRp)
    warning('co2:NoRpCase', 'Late-core case %s / %s not found; nominal r_p only.', ...
            cfg.rp.case.Mixture, cfg.rp.case.Rock);
else
    cRp = dp.getCase(iRp(1));
    fprintf('=== Stage B: late-core r_p from %s %s (shift %.1f +/- %.1f psi) ===\n', ...
            cRp.Mixture, cRp.Rock, dp.ShiftMid_psi(iRp(1)), dp.ConfHalfWidth_psi(iRp(1)));
    for f = frozen
        runner.setFWI(f.k, f.p, f.lambda);
        runner.clearWarm();
        obj = @(rp_nm) rpResidual(rp_nm, runner, dp, iRp(1), cfg, f, cfg.rp.mode);
        rpOpt = fminbnd(obj, cfg.rp.bracket_nm(1), cfg.rp.bracket_nm(2), optimset('TolX', 1e-3));
        [resid, shift, tip] = rpResidual(rpOpt, runner, dp, iRp(1), cfg, f, cfg.rp.mode);
        rpDerived(char(f.name)) = rpOpt * 1e-9;
        rpRows(end+1, :) = {f.name, f.role, rpOpt, resid, shift, dp.ShiftMid_psi(iRp(1)), ...
            dp.ConfHalfWidth_psi(iRp(1)), tip, cfg.rp.band_nm(1), cfg.rp.band_nm(2)}; %#ok<AGROW>
        fprintf('  %-20s derived r_p = %.2f nm (shift %.1f vs %.1f psi) | band carried [%.2f, %.2f] nm\n', ...
                f.name, rpOpt, shift, dp.ShiftMid_psi(iRp(1)), cfg.rp.band_nm);
    end
    rpTable = cell2table(rpRows, 'VariableNames', {'Model', 'Role', 'rp_derived_nm', 'Residual_psi', ...
        'Shift_model_psi', 'Shift_exp_psi', 'HalfWidth_psi', 'TIP_used', 'Band_lo_nm', 'Band_hi_nm'});
    rpTable = convertvars(rpTable, {'Model', 'Role'}, 'string');
    fprintf('\n');
end

%% 4. Stage C: bulk CO2 dew points with literature BIPs ------------------------
fprintf('=== Stage C: bulk CO2 dew points, literature BIPs, residual reported as found ===\n');
bulkRows = {};
for i = find(isCO2).'
    c = dp.getCase(i);
    [Pb, ~, okb] = runner.bulk(c);
    bulkRows(end+1, :) = {c.Mixture, c.Rock, strjoin(string(c.z.'), ","), dp.BulkPdew_psi(i), ...
        Pb / psi, Pb / psi - dp.BulkPdew_psi(i), 100 * (Pb / (dp.BulkPdew_psi(i) * psi) - 1), okb}; %#ok<AGROW>
end
bulkTable = cell2table(bulkRows, 'VariableNames', {'Mixture', 'Rock', 'z', 'Bulk_exp_psi', ...
    'Bulk_model_psi', 'Residual_psi', 'Residual_pct', 'Converged'});
bulkTable = convertvars(bulkTable, {'Mixture', 'Rock', 'z'}, 'string');
disp(bulkTable);
fprintf(['The BIPs are the literature values in the BIP sheet and the compositions follow the ' ...
         'stated policy; the residuals above are reported as found, not tuned away.\n\n']);

%% 5. Stages D and E: confined CO2 predictions ---------------------------------
predRows = {};  pairRows = {};
for f = frozen
    runner.setFWI(f.k, f.p, f.lambda);
    for capMode = [false, true]
        if capMode && ~cfg.tip.capAtCalibrationEdge, continue; end
        tipLabel = "raw extrapolation";
        if capMode, tipLabel = "M capped at calibration edge"; end

        for rpLabel = ["nominal", "derived"]
            for md = cfg.modes
                runner.clearWarm();
                for i = find(isCO2 & dp.HasConfined).'
                    c = dp.getCase(i);
                    rp = poreRadiusFor(c, rpLabel, cfg, rpDerived, f);
                    if ~isfinite(rp), continue; end
                    eng = runner.engine(c.Mixture, c.Rock);
                    [K, pairs] = tipMatrix(eng, c, rp, f, cfg, capMode);
                    for q = 1:height(pairs)
                        pairRows(end+1, :) = {f.name, tipLabel, rpLabel, c.Mixture, c.Rock, ...
                            rp * 1e9, pairs.Pair(q), pairs.M(q), pairs.M_used(q), ...
                            pairs.sigma_over_rp(q), pairs.TIP(q)}; %#ok<AGROW>
                    end
                    if any(K(:) >= 0 & K(:) ~= 0)
                        error('co2:PositiveTIP', 'Correlation produced TIP >= 0 for %s.', c.Mixture);
                    end
                    r = runner.solve(c, md, K, 'r', rp, ...
                            'Key', sprintf("%s|%d|%s|%s|%d", f.name, i, md, rpLabel, capMode));
                    err = r.shift_psi - dp.ShiftMid_psi(i);
                    predRows(end+1, :) = {f.name, f.role, md, tipLabel, rpLabel, rp * 1e9, ...
                        c.Mixture, c.Rock, dp.ShiftMid_psi(i), dp.ConfHalfWidth_psi(i), ...
                        r.shift_psi, err, abs(err) <= dp.ConfHalfWidth_psi(i), ...
                        100 * (r.P / (dp.ConfMid_psi(i) * psi) - 1), ...
                        abs(100 * (r.P / (dp.ConfMid_psi(i) * psi) - 1)) <= cfg.accept.maxPdewErr_pct, ...
                        r.ok, r.route, r.reason}; %#ok<AGROW>
                end
            end
        end
    end
end
predTable = cell2table(predRows, 'VariableNames', {'Model', 'Role', 'Mode', 'TIPsetting', 'rpSetting', ...
    'rp_nm', 'Mixture', 'Rock', 'Shift_exp_psi', 'HalfWidth_psi', 'Shift_pred_psi', 'Error_psi', ...
    'InBand', 'Pdew_err_pct', 'Within1pct', 'Accepted', 'Route', 'Reason'});
pairTable = cell2table(pairRows, 'VariableNames', {'Model', 'TIPsetting', 'rpSetting', 'Mixture', ...
    'Rock', 'rp_nm', 'Pair', 'M_geomean', 'M_used', 'sigma_over_rp', 'TIP'});
for nm = ["predTable", "pairTable"]
    t = eval(nm);
    isTxt = varfun(@(v) iscell(v) && all(cellfun(@(e) isstring(e) || ischar(e), v)), t, ...
                   'OutputFormat', 'uniform');
    if any(isTxt), t = convertvars(t, t.Properties.VariableNames(isTxt), 'string'); end
    eval(nm + " = t;");
end

fprintf('=== Stages D/E: confined CO2 predictions ===\n');   disp(predTable);
fprintf('=== Pair TIPs used (raw vs capped) ===\n');          disp(unique(pairTable, 'stable'));
co2c1 = pairTable(contains(pairTable.Pair, "CO2") & contains(pairTable.Pair, "C1") & ...
                  ~contains(pairTable.Pair, "nC"), :);
if ~isempty(co2c1)
    fprintf(['\nCO2-C1 sensitivity: M_geomean = %.2f lies below the calibration range ' ...
             '(cap %.2f). Raw and capped TIPs are both reported above.\n'], ...
            co2c1.M_geomean(1), cfg.tip.Mcap);
end

%% 6. Freeze, tag and write ------------------------------------------------------
frozenTable = struct2table(frozen);
frozenTable.Tag = repmat(cfg.freezeTag, height(frozenTable), 1);
frozenTable.TIPform = repmat("TIP = -exp(alpha + beta/sqrt(Mi*Mj) + gamma*sigma_bar/rp)", ...
                             height(frozenTable), 1);
frozenTable.Normalisation = repmat("blocks scaled by RMS of their own data; w on the dTc block", ...
                             height(frozenTable), 1);
frozenTable.Mdef = repmat("geomean", height(frozenTable), 1);
frozenTable.w = [NaN; 1];
frozenTable.CompositionPolicy = repmat(cfg.composition.policy, height(frozenTable), 1);

out = fullfile(cfg.resultsDir, cfg.outFile);
save(out + ".mat", 'cfg', 'frozen', 'frozenTable', 'bulkTable', 'predTable', 'pairTable');
if exist('tcTable', 'var'), save(out + ".mat", 'tcTable', '-append'); end
if exist('rpTable', 'var'), save(out + ".mat", 'rpTable', '-append'); end
writetable(frozenTable, fullfile(cfg.resultsDir, "frozen_model.csv"));
xls = out + ".xlsx";
if isfile(xls), delete(xls); end
writetable(frozenTable, xls, 'Sheet', 'Frozen_model');
if exist('tcTable', 'var'), writetable(tcTable, xls, 'Sheet', 'A_CO2_Tc_blind'); end
if exist('rpTable', 'var'), writetable(rpTable, xls, 'Sheet', 'B_rp_derivation'); end
writetable(bulkTable, xls, 'Sheet', 'C_bulk_CO2');
writetable(predTable, xls, 'Sheet', 'D_predictions');
writetable(pairTable, xls, 'Sheet', 'E_pair_TIPs');
fprintf('\nWritten: %s\n         %s\n', xls, fullfile(cfg.resultsDir, "frozen_model.csv"));

%% ============================== Local functions ==============================
function rp = poreRadiusFor(c, rpLabel, cfg, rpDerived, f)
    % Nominal: the sheet's own r_p column if present, else the rock default.
    % Derived: the late-core r_p from stage B for this frozen model.
    if rpLabel == "derived"
        rp = NaN;
        if isKey(rpDerived, char(f.name)), rp = rpDerived(char(f.name)); end
        return;
    end
    rp = c.PoreRadius_m;
    if ~isfinite(rp), rp = cfg.poreRadiusNominal(char(c.Rock)); end
end

function [K, T] = tipMatrix(eng, c, rp, f, cfg, capM)
    % Every pair TIP from the frozen correlation, with optional capping of M at
    % the calibration edge (declared sensitivity for CO2-C1).
    sig = eng.eos.Fluid.LJ_Size(:);
    MW  = eng.eos.Fluid.MW(:);
    nm  = eng.eos.Fluid.ComponentNames;
    n   = numel(c.z);
    M   = sqrt(MW(:) * MW(:).');
    Mu  = M;
    if capM, Mu = max(M, cfg.tip.Mcap); end
    sR  = 0.5 * (sig(:) + sig(:).') / rp;
    K   = -exp(f.alpha + f.beta ./ Mu + f.gamma * sR);
    K(1:n+1:end) = 0;
    K   = (K + K.') / 2;
    rows = {};
    for a = 1:n-1
        for b = a+1:n
            rows(end+1, :) = {nm(a) + "-" + nm(b), M(a, b), Mu(a, b), sR(a, b), K(a, b)}; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'Pair', 'M', 'M_used', 'sigma_over_rp', 'TIP'});
    T.Pair = string(T.Pair);
end

function [resid, shift, tip] = rpResidual(rp_nm, runner, dp, i, cfg, f, mode)
    % |model shift - experimental shift| at pore radius rp_nm, with the TIP taken
    % from the frozen correlation at that same r_p (self-consistent).
    c   = dp.getCase(i);
    rp  = rp_nm * 1e-9;
    eng = runner.engine(c.Mixture, c.Rock);
    [K, T] = tipMatrix(eng, c, rp, f, cfg, false);
    tip = T.TIP(1);
    r = runner.solve(c, mode, K, 'r', rp, 'Key', sprintf("rp|%s|%.4f", f.name, rp_nm));
    shift = NaN;  resid = 1e6;
    if r.ok
        shift = r.shift_psi;
        resid = abs(shift - dp.ShiftMid_psi(i));
    end
end