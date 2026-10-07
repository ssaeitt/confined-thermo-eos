%% CO2_PIPELINE_V3  CO2-HC mixtures with the frozen model, rock-offset TIP (FWI-Only)
%
% v3 changes relative to v2
%   * r_p on the SHIFT basis (model confined - model bulk), cfg.rpSource = "manual" by default.
%     v2 read the absolute-basis r_p from backcalc_rp_B1_late.mat (anchor shift 239.6 psi,
%     outside its 216..235 band).
%   * Stage E: exact-match BIP fronts for the unreconciled bulk dew points. Along each front the
%     CO2-HC BIPs are stepped and, at every step, each CO2 mixture's CO2 fraction (C1:nC8 held)
%     is solved so that its model bulk Pdew equals the measured bulk Pdew exactly. The confined
%     shift is then predicted at that (BIP, composition) pair.
%        front "CO2C1": BIP CO2-C1 stepped, BIP CO2-nC8 fixed at the literature value
%        front "both" : BIP CO2-C1 and BIP CO2-nC8 raised by the same increment
%     BIPs are shared by both mixtures (same pairs); compositions are per mixture. Each front point
%     reports the share of the nominal-composition bulk mismatch absorbed by the BIP change
%     (BIPShare: 0 = composition only, 1 = BIP only). The BIP-only end of each front (nominal
%     composition) is located exactly per mixture.
%     BIPs are changed by writing a temporary copy of the mixture workbook (BIP sheet) and building
%     a fresh DewPointRunner from it, so no solver internals are touched. C1-nC8 is unaffected,
%     so the r_p anchor and r_p are unchanged.
%   Terminology: "BIP" = binary interaction parameter; "k" is only the FWI parameter.
%
% Supersedes CO2_pipeline_v2.m. Nothing in the confinement model is fitted here.
%
% Frozen model
%   k, p, lambda       : per model (jL primary; w = 1 optional)
%   TIP_ij             : -exp( alpha_rock + B * x_ij )      (sheet "Correlation" of the TIP workbook)
%   primary            : jL / LL-inv_wbar   x = 1/wbar_ij, wbar = (omega_i + omega_j)/2
%   sensitivity        : jL / LL-inv_Tstar  x = sqrt(eps_i eps_j)/(kB T)   (no omega)
%   optional           : pareto-exact-w1 / LL-Tstar  (cfg.models(3).enabled; needs its own r_p, see below)
%   r_p (late B1 core) : read from results/backcalc_rp_B1_late.mat (Backcalc_rp_B1_late.m),
%                        manual fallback mid 7.82 nm, band edges 7.61 / 8.05 nm
%
% Stages
%   A  blind CO2 critical-shift test (depends on k, p, lambda only)
%   B  base case: nominal composition, literature BIPs, r_p = derived mid value; bulk residual
%      and shift reported. The 0 % CO2 row is the r_p anchor (matched by construction); the
%      20 % and 40 % mixtures are the genuine tests.
%   C  one-at-a-time variations of the base case:
%        effcomp     effective CO2 fraction that reproduces the experimental bulk Pdew
%        rp_lo/rp_hi r_p at the band edges
%        rp_nominal  r_p = 11.25 nm (same-core hypothesis)
%        xclip       descriptor x clipped to the calibrated range (replaces the old M cap)
%   D  success criteria (unchanged): both shifts > 0; shift(40 %) < shift(20 %); both shifts near
%      their measured bands (20 %: +228..+247 psi, 40 %: +141..+159 psi).
% Shift = model confined Pdew - model bulk Pdew (headline). Absolute Pdew error is non-gating.

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

%% 0. Configuration ---------------------------------------------------------------
cfg.configDir   = fullfile(scriptDir, "config");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.mixtureFile = fullfile(cfg.configDir, "MixtureData - Complete Mixtures.xlsx");
cfg.tipFile     = fullfile(cfg.configDir, "tip_correlation_rockoffset.xlsx");
cfg.tcShiftFile = fullfile(cfg.configDir, "Extended_Dataset - Original.xlsx");
cfg.tcSheets    = ["All Data (All comp)", "All Data"];     % first sheet with CO2 rows is used
cfg.outFile     = "co2_pipeline_v3";
cfg.T           = 293.15;
cfg.mode        = "FWIOnly";
cfg.kB          = 1.380649e-23;
cfg.wallRock    = "EF2";                % row used only to read mineral_E (Tc calibration wall)
cfg.wallMineral = "Quartz";
cfg.poreRadiusNominal = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
cfg.run.tc = true;                      % stage A can run on its own (set the rest off to test it)
cfg.run.mixtures = true;

% Models. rpFile = Backcalc_rp_B1_late.m output; rpManual = [mid lo hi] in nm used if the file is absent.
% w = 1 is OFF by default: LL-inv_wbar failed the criteria with w = 1, and the Tstar forms that passed
% need their own r_p (run Backcalc_rp_B1_late.m with cfg.preset = "w1", then set enabled = true).
cfg.models = struct( ...
    'name',        {"jL-inv_wbar", "jL-inv_Tstar", "w1-Tstar"}, ...
    'role',        {"primary", "sensitivity", "optional"}, ...
    'combination', {"joint-linearized", "joint-linearized", "pareto-exact-w1"}, ...
    'form',        {"LL-inv_wbar", "LL-inv_Tstar", "LL-Tstar"}, ...
    'k',           {112.42634826677, 112.42634826677, 166.410385934507}, ...
    'p',           {1.3558179290678, 1.3558179290678, 1.5488421885245}, ...
    'lambda',      {2.0916235090532, 2.0916235090532, 4.47908205631447}, ...
    'rpFile',      {"backcalc_rp_B1_late.mat", "backcalc_rp_B1_late.mat", "backcalc_rp_B1_late_w1.mat"}, ...
    'rpManual',    {[7.82711826617637 7.62188369515433 8.04796163088063], [7.82711826617637 7.62188369515433 8.04796163088063], [8.7824 8.6 8.982]}, ...
    'enabled',     {true, true, true});
cfg.tipMode = "FWIOnly";
% r_p source: "manual" uses rpManual = [mid lo hi] nm on the SHIFT basis (project convention:
% shift = model confined - model bulk). Both jL forms share one r_p: the anchor pair C1-nC8 on B1
% has the same TIP in every form. w1 has its own r_p (different k, p, lambda and anchor TIP).
% "file" reads results/<rpFile> (v2 behaviour; that file held the absolute-basis value).
cfg.rpSource = "manual";

% Stage E: exact-match BIP fronts (BIP = binary interaction parameter of the PR a_m)
cfg.front.run      = true;
cfg.front.bipLit   = struct('CO2_C1', 0.09897, 'CO2_nC8', 0.11541);   % literature (BIP sheet)
cfg.front.CO2C1    = [0.09897, 0.15:0.05:0.60];   % BIP CO2-C1 values; BIP CO2-nC8 fixed
cfg.front.both     = 0:0.05:0.50;                 % common increment added to both CO2 BIPs
cfg.front.rp       = "mid";
cfg.front.clip     = false;                       % TIP as in the base case
cfg.front.workDir  = fullfile(scriptDir, "results", "bip_front_tmp");
cfg.front.checkTol_psi = 0.5;                     % literature point must reproduce the base bulk

% Scenarios (one change at a time relative to "base")
cfg.scen = struct( ...
    'name', {"base", "effcomp", "rp_lo", "rp_hi", "rp_nominal", "xclip"}, ...
    'comp', {"nominal", "effective", "nominal", "nominal", "nominal", "nominal"}, ...
    'rp',   {"mid", "mid", "lo", "hi", "nominal", "mid"}, ...
    'clip', {false, false, false, false, false, true});

% Anchor row (sets r_p) and criteria
cfg.anchor = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.criteria.nearMargin_psi = 25;       % "near the band" = within the band widened by this margin (assumption)
cfg.composition.policy = "BIPs are fixed at literature values and compositions are used as charged (nominal). " + ...
    "The effective-composition scenario is a diagnostic only; unreconciled bulk dew points are reported as such.";
cfg.composition.scanGrid = 0:0.02:0.6;  % CO2 mole fractions scanned for the effective-CO2 inversion
cfg.composition.xRound   = 4;           % runner bulk cache keys keep 5 significant digits

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
models = cfg.models([cfg.models.enabled]);
psi = entities.DewPointDataset.psi2Pa;
fprintf('Models: %s\n%s\n\n', strjoin([models.name], ", "), cfg.composition.policy);

%% 1. Data and engines -------------------------------------------------------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
isCO2 = contains(dp.Mixture, "CO2") & dp.HasConfined;
iAnch = dp.findCases(cfg.anchor.Mixture, cfg.anchor.Rock, cfg.anchor.z);
if numel(iAnch) ~= 1
    error('co2:AnchorNotUnique', '%d rows match the r_p anchor; expected one.', numel(iAnch));
end
iCO2 = find(isCO2).';
if numel(iCO2) < 2
    error('co2:NeedTwoCO2', 'Need the 20 %% and 40 %% CO2 rows with confined data (found %d).', numel(iCO2));
end
xco2 = zeros(size(iCO2));
for q = 1:numel(iCO2)
    cq = dp.getCase(iCO2(q));  xco2(q) = cq.z(find(cq.Components == "CO2", 1));
end
[xco2, ord] = sort(xco2);  iCO2 = iCO2(ord);
iLow = iCO2(1);  iHigh = iCO2(end);                  % 20 % and 40 % rows
caseIdx = [iAnch, iCO2];

runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadiusNominal, 'maxShift_psi', Inf);
for i = caseIdx, runner.engine(dp.Mixture(i), dp.Rock(i)); end

%% 2. Stage A: blind CO2 critical-shift test ---------------------------------------
tcTable = table();
if cfg.run.tc
    fluidAll = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile);
    rockRef  = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallRock);
    q = rockRef.mineral(cfg.wallMineral);
    calRock = entities.RockProperties("CalWall-" + q.Name, 0, q.Name, q.Epsilon_K, q.GrainDensity, 1);
    eosCal = thermo.ConfinedEOS(fluidAll, calRock);

    tcCO2 = [];  sheetUsed = "";
    for sh = cfg.tcSheets
        try
            Tt = entities.TcShiftDataset.loadFromWorkbook(cfg.tcShiftFile, 'Sheet', sh, 'Components', "CO2");
        catch
            continue;
        end
        if Tt.NumPoints > 0, tcCO2 = Tt;  sheetUsed = sh;  break; end
    end
    if isempty(tcCO2)
        warning('co2:NoCO2Tc', 'No CO2 critical-shift points found in %s. Stage A skipped.', cfg.tcShiftFile);
    else
        [idx, r] = tcCO2.modelInputs(fluidAll);
        kpl = unique([[models.k].', [models.p].', [models.lambda].'], 'rows', 'stable');
        fprintf('=== Stage A: blind CO2 critical-shift test (%d points, sheet "%s") ===\n', tcCO2.NumPoints, sheetUsed);
        rowsA = {};
        for u = 1:size(kpl, 1)
            lab = strjoin([models(ismember([[models.k].', [models.p].', [models.lambda].'], kpl(u, :), 'rows')).name], "+");
            eosCal.k = kpl(u, 1);  eosCal.pT_wall = kpl(u, 2);  eosCal.lambda = kpl(u, 3);
            for n = 1:tcCO2.NumPoints
                i = idx(n);  pred = NaN;
                try
                    pred = 1 - eosCal.pureCriticalPoint(i, r(n)) / eosCal.Fluid.Tc(i);
                catch ME
                    if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
                end
                rowsA(end+1, :) = {lab, tcCO2.Component(n), tcCO2.dH(n), tcCO2.Reference(n), tcCO2.dTc(n), ...
                    pred, pred - tcCO2.dTc(n), 100 * (pred - tcCO2.dTc(n)) / tcCO2.dTc(n)}; %#ok<AGROW>
            end
        end
        tcTable = textToString(cell2table(rowsA, 'VariableNames', {'Params', 'Component', 'dH', 'Reference', ...
            'dTc_exp', 'dTc_pred', 'Residual', 'Error_pct'}));
        disp(tcTable);
        fprintf('Blind test: CO2 was excluded from the HC-only calibration data (depends on k, p, lambda only).\n\n');
    end
end
if ~cfg.run.mixtures
    fprintf('cfg.run.mixtures = false: stopping after stage A.\n');
    return;
end

%% 3. TIP coefficients, calibrated x range, r_p per model ---------------------------
tipTab  = readtable(cfg.tipFile, 'Sheet', 'Correlation', 'VariableNamingRule', 'preserve');
pairTab = readtable(cfg.tipFile, 'Sheet', 'Pair_TIPs',   'VariableNamingRule', 'preserve');
for m = 1:numel(models)
    mm = models(m);
    sel = string(pairTab.Combination) == mm.combination & string(pairTab.Mode) == cfg.tipMode & ...
          string(pairTab.Form) == mm.form & string(pairTab.Context) == "binary";
    if ~any(sel), error('co2:NoPairRows', 'No binary Pair_TIPs rows for %s / %s.', mm.combination, mm.form); end
    models(m).xRange = [min(pairTab.x(sel)), max(pairTab.x(sel))];
    models(m).rp = loadRp(mm, cfg);
    fprintf('%-13s (%s): %s / %s | calibrated x range [%.4f, %.4f] | r_p mid %.3f nm, edges [%.3f, %.3f] (%s)\n', ...
        mm.name, mm.role, mm.combination, mm.form, models(m).xRange, models(m).rp.mid, models(m).rp.lo, ...
        models(m).rp.hi, models(m).rp.src);
    if ~isfinite(models(m).rp.mid)
        error('co2:NoRp', 'No r_p for %s. Run Backcalc_rp_B1_late.m (preset for this model) first.', mm.name);
    end
end
fprintf('\n');

%% 4. Experimental references, model bulk, effective CO2 ----------------------------
nCase = numel(caseIdx);
ref = struct('bulkExp', num2cell(dp.BulkPdew_psi(caseIdx).'), 'shLo', num2cell((dp.ConfLower_psi(caseIdx) - dp.BulkPdew_psi(caseIdx)).'), ...
             'shUp', num2cell((dp.ConfUpper_psi(caseIdx) - dp.BulkPdew_psi(caseIdx)).'), ...
             'confLo', num2cell(dp.ConfLower_psi(caseIdx).'), 'confUp', num2cell(dp.ConfUpper_psi(caseIdx).'));
bulkModel = nan(1, nCase);  zEff = cell(1, nCase);  xEff = nan(1, nCase);
fprintf('=== Bulk dew points (literature BIPs, charged composition) ===\n');
for a = 1:nCase
    c = dp.getCase(caseIdx(a));
    [Pb, ~, ok] = runner.bulk(c);
    if ok, bulkModel(a) = Pb / psi; end
    resid = 100 * (bulkModel(a) / ref(a).bulkExp - 1);
    msg = "";
    if caseIdx(a) ~= iAnch
        [xr, ~, ~] = effectiveCO2(runner, c, ref(a).bulkExp, cfg);
        if ~isempty(xr)
            iC = find(c.Components == "CO2", 1);
            [~, j] = min(abs(xr - c.z(iC)));
            xEff(a) = xr(j);
            rest = c.z;  rest(iC) = 0;  rest = rest / sum(rest);
            z = (1 - xEff(a)) * rest;  z(iC) = xEff(a);  zEff{a} = z(:);
            msg = sprintf(' | effective CO2 %.2f mol%% (nominal %.1f)', 100 * xEff(a), 100 * c.z(iC));
        else
            msg = ' | no CO2 fraction reproduces the experimental bulk';
        end
    end
    fprintf('  %-12s z = [%s]: exp %.0f, model %.1f psi (%+.2f %%)%s\n', c.Mixture, num2str(c.z.', '%.3f '), ...
        ref(a).bulkExp, bulkModel(a), resid, msg);
end
fprintf('\n');

%% 5. Scenarios ----------------------------------------------------------------------
runRows = {};  tipRows = {};
for m = 1:numel(models)
    mm = models(m);
    runner.setFWI(mm.k, mm.p, mm.lambda);
    runner.clearWarm();
    for sc = cfg.scen
        for a = 1:nCase
            i = caseIdx(a);
            c = dp.getCase(i);
            isAnch = (i == iAnch);
            if sc.comp == "effective"
                if isAnch || isempty(zEff{a}), continue; end
                c.z = zEff{a};
            end
            switch sc.rp
                case "mid",     rpnm = mm.rp.mid;
                case "lo",      rpnm = mm.rp.lo;
                case "hi",      rpnm = mm.rp.hi;
                case "nominal", rpnm = 1e9 * cfg.poreRadiusNominal(char(c.Rock));
            end
            eng = runner.engine(c.Mixture, c.Rock);
            co  = tipCoefficients(tipTab, mm.combination, cfg.tipMode, mm.form, c.Rock);
            [K, pt] = buildTIP(eng.eos.Fluid, co, cfg.T, cfg.kB, mm.xRange, sc.clip);
            for q = 1:height(pt)
                tipRows(end+1, :) = {mm.name, sc.name, c.Mixture, strjoin(compose("%.4f", c.z(:).'), ","), ...
                    pt.Pair(q), pt.x(q), pt.x_used(q), pt.RangeFlag(q), pt.TIP(q)}; %#ok<AGROW>
            end
            r = runner.solve(c, cfg.mode, K, 'r', rpnm * 1e-9, ...
                'Key', mm.name + "|" + sc.name + "|" + i);
            Pb_scn = r.PdewBulkModel / psi;                     % model bulk at the scenario composition
            shift = r.shift_psi;
            [inBand, miss] = bandMiss(shift, ref(a).shLo, ref(a).shUp);
            pdErr = 100 * (r.P / psi / (0.5 * (ref(a).confLo + ref(a).confUp)) - 1);
            runRows(end+1, :) = {mm.name, mm.role, sc.name, c.Mixture, strjoin(compose("%.4f", c.z(:).'), ","), ...
                isAnch, rpnm, ref(a).bulkExp, Pb_scn, 100 * (Pb_scn / ref(a).bulkExp - 1), r.P / psi, shift, ...
                ref(a).shLo, ref(a).shUp, r.shiftExpBulk_psi, inBand, miss, pdErr, r.ok, string(r.route), string(r.reason)}; %#ok<AGROW>
        end
    end
end
runTable = textToString(cell2table(runRows, 'VariableNames', {'Model', 'Role', 'Scenario', 'Mixture', 'z', ...
    'IsRpAnchor', 'rp_nm', 'Bulk_exp_psi', 'Bulk_model_psi', 'BulkResid_pct', 'Pconf_model_psi', ...
    'Shift_modelBulk_psi', 'ShiftExp_lo_psi', 'ShiftExp_up_psi', 'Shift_vsExpBulk_psi', 'InBand', ...
    'MissOutsideBand_psi', 'PdewErr_pct_nonGating', 'Accepted', 'Route', 'Reason'}));
tipTable = textToString(cell2table(tipRows, 'VariableNames', {'Model', 'Scenario', 'Mixture', 'z', 'Pair', ...
    'x', 'x_used', 'RangeFlag', 'TIP'}));
tipTable = unique(tipTable, 'stable');

%% 6. Success criteria ------------------------------------------------------------------
critRows = {};
for m = 1:numel(models)
    for sc = cfg.scen
        selM = runTable.Model == models(m).name & runTable.Scenario == sc.name & ~runTable.IsRpAnchor;
        T0 = runTable(selM, :);
        if height(T0) < 2, continue; end
        [~, o] = sort(extractXco2(T0.z));
        sLow = T0.Shift_modelBulk_psi(o(1));  sHigh = T0.Shift_modelBulk_psi(o(end));
        loL = T0.ShiftExp_lo_psi(o(1));   upL = T0.ShiftExp_up_psi(o(1));
        loH = T0.ShiftExp_lo_psi(o(end)); upH = T0.ShiftExp_up_psi(o(end));
        [~, missL] = bandMiss(sLow, loL, upL);  [~, missH] = bandMiss(sHigh, loH, upH);
        c1 = sLow > 0 && sHigh > 0;
        c2 = sHigh < sLow;
        c3 = abs(missL) <= cfg.criteria.nearMargin_psi && abs(missH) <= cfg.criteria.nearMargin_psi;
        critRows(end+1, :) = {models(m).name, models(m).role, sc.name, sLow, sHigh, missL, missH, c1, c2, c3, c1 && c2 && c3}; %#ok<AGROW>
    end
end
critTable = textToString(cell2table(critRows, 'VariableNames', {'Model', 'Role', 'Scenario', 'Shift_low_psi', ...
    'Shift_high_psi', 'MissLow_psi', 'MissHigh_psi', 'BothPositive', 'HighLessThanLow', 'NearBands', 'AllPass'}));

fprintf('=== Anchor (0 %% CO2, matched by construction) ===\n');
disp(runTable(runTable.IsRpAnchor & ismember(runTable.Scenario, ["base", "rp_lo", "rp_hi", "rp_nominal"]), ...
    {'Model', 'Scenario', 'rp_nm', 'Pconf_model_psi', 'Shift_modelBulk_psi', 'ShiftExp_lo_psi', 'ShiftExp_up_psi', 'InBand'}));
fprintf('=== CO2 mixtures: shifts vs measured bands (20 %%: +228..+247, 40 %%: +141..+159 psi) ===\n');
disp(runTable(~runTable.IsRpAnchor, {'Model', 'Scenario', 'Mixture', 'z', 'rp_nm', 'BulkResid_pct', ...
    'Shift_modelBulk_psi', 'ShiftExp_lo_psi', 'ShiftExp_up_psi', 'InBand', 'MissOutsideBand_psi', 'Accepted'}));
fprintf('=== Success criteria (near = band +/- %g psi) ===\n', cfg.criteria.nearMargin_psi);
disp(critTable);
fprintf('=== Pair TIPs used (raw and clipped) ===\n');
disp(tipTable(tipTable.Mixture ~= cfg.anchor.Mixture, :));

%% 7. Stage E: exact-match BIP fronts ---------------------------------------------------
frontPts = table();  frontRuns = table();  frontCrit = table();  bipOnly = table();
if cfg.front.run
    if ~isfolder(cfg.front.workDir), mkdir(cfg.front.workDir); end
    iCO2pos = find(caseIdx ~= iAnch);                       % positions in caseIdx of the CO2 rows
    fronts = struct('name', {"CO2C1", "both"}, 'grid', {cfg.front.CO2C1, cfg.front.both});
    ptRows = {};  runRowsE = {};
    PnomGrid = struct();
    for fr = fronts
        nG = numel(fr.grid);
        PnomGrid.(fr.name) = nan(numel(iCO2pos), nG);
        bipC1 = nan(1, nG);  bipC8 = nan(1, nG);
        for g = 1:nG
            [bipC1(g), bipC8(g)] = frontBIPs(fr.name, fr.grid(g), cfg);
            bipFile = fullfile(cfg.front.workDir, sprintf("bip_%s_%02d.xlsx", fr.name, g));
            writeBipWorkbook(cfg.mixtureFile, bipFile, {"CO2", "C1", bipC1(g); "CO2", "nC8", bipC8(g)});
            rr = solvers.DewPointRunner(bipFile, cfg.poreRadiusNominal, 'maxShift_psi', Inf);
            for i = caseIdx, rr.engine(dp.Mixture(i), dp.Rock(i)); end
            fprintf('Front %s, point %d/%d: BIP CO2-C1 = %.4f, BIP CO2-nC8 = %.4f\n', fr.name, g, nG, bipC1(g), bipC8(g));
            for qa = 1:numel(iCO2pos)
                a = iCO2pos(qa);  c = dp.getCase(caseIdx(a));
                iC = find(c.Components == "CO2", 1);
                [Pn, ~, okn] = rr.bulk(c);
                Pnom = NaN;  if okn, Pnom = Pn / psi; end
                PnomGrid.(fr.name)(qa, g) = Pnom;
                if abs(bipC1(g) - cfg.front.bipLit.CO2_C1) < 1e-12 && abs(bipC8(g) - cfg.front.bipLit.CO2_nC8) < 1e-12 ...
                        && ~(abs(Pnom - bulkModel(a)) <= cfg.front.checkTol_psi)
                    error('co2:BipFileCheck', ['The temporary-workbook runner does not reproduce the base bulk ' ...
                        'at literature BIPs (%s: %.2f vs %.2f psi).'], c.Mixture, Pnom, bulkModel(a));
                end
                share = (Pnom - bulkModel(a)) / (ref(a).bulkExp - bulkModel(a));
                xr = effectiveCO2(rr, c, ref(a).bulkExp, cfg);
                xE = NaN;  zE = [];
                if ~isempty(xr)
                    [~, j] = min(abs(xr - c.z(iC)));  xE = xr(j);
                    rest = c.z;  rest(iC) = 0;  rest = rest / sum(rest);
                    zE = (1 - xE) * rest;  zE(iC) = xE;  zE = zE(:);
                end
                ptRows(end+1, :) = {fr.name, g, bipC1(g), bipC8(g), c.Mixture, c.z(iC), ref(a).bulkExp, ...
                    Pnom, share, xE, 100 * (xE - c.z(iC)), numel(xr)}; %#ok<AGROW>
                if isempty(zE), continue; end
                cE = c;  cE.z = zE;
                for m = 1:numel(models)
                    mm = models(m);
                    rr.setFWI(mm.k, mm.p, mm.lambda);
                    rpnm = mm.rp.(cfg.front.rp);
                    eng = rr.engine(cE.Mixture, cE.Rock);
                    co  = tipCoefficients(tipTab, mm.combination, cfg.tipMode, mm.form, cE.Rock);
                    K   = buildTIP(eng.eos.Fluid, co, cfg.T, cfg.kB, mm.xRange, cfg.front.clip);
                    r = rr.solve(cE, cfg.mode, K, 'r', rpnm * 1e-9, ...
                        'Key', "front|" + fr.name + "|" + g + "|" + mm.name + "|" + caseIdx(a));
                    sh = NaN;  if r.ok, sh = r.shift_psi; end
                    [inB, miss] = bandMiss(sh, ref(a).shLo, ref(a).shUp);
                    runRowsE(end+1, :) = {fr.name, g, bipC1(g), bipC8(g), mm.name, mm.role, c.Mixture, ...
                        c.z(iC), xE, share, rpnm, r.PdewBulkModel / psi, r.P / psi, sh, ref(a).shLo, ref(a).shUp, ...
                        inB, miss, r.ok, string(r.route), string(r.reason)}; %#ok<AGROW>
                end
            end
        end
    end
    frontPts = textToString(cell2table(ptRows, 'VariableNames', {'Front', 'Point', 'BIP_CO2_C1', 'BIP_CO2_nC8', ...
        'Mixture', 'xCO2_nominal', 'Bulk_exp_psi', 'Bulk_model_nominalComp_psi', 'BIPShare', ...
        'xCO2_exactMatch', 'dCO2_molpct', 'NumRoots'}));
    frontRuns = textToString(cell2table(runRowsE, 'VariableNames', {'Front', 'Point', 'BIP_CO2_C1', 'BIP_CO2_nC8', ...
        'Model', 'Role', 'Mixture', 'xCO2_nominal', 'xCO2_exactMatch', 'BIPShare', 'rp_nm', 'Bulk_model_psi', ...
        'Pconf_model_psi', 'Shift_modelBulk_psi', 'ShiftExp_lo_psi', 'ShiftExp_up_psi', 'InBand', ...
        'MissOutsideBand_psi', 'Accepted', 'Route', 'Reason'}));

    % criteria per front point and model (same definitions as stage D)
    cr = {};
    for fr = fronts
        for g = 1:numel(fr.grid)
            for m = 1:numel(models)
                T0 = frontRuns(frontRuns.Front == fr.name & frontRuns.Point == g & frontRuns.Model == models(m).name, :);
                if height(T0) < 2, continue; end
                [~, o] = sort(T0.xCO2_nominal);
                sL = T0.Shift_modelBulk_psi(o(1));  sH = T0.Shift_modelBulk_psi(o(end));
                [~, mL] = bandMiss(sL, T0.ShiftExp_lo_psi(o(1)), T0.ShiftExp_up_psi(o(1)));
                [~, mH] = bandMiss(sH, T0.ShiftExp_lo_psi(o(end)), T0.ShiftExp_up_psi(o(end)));
                c1 = sL > 0 && sH > 0;  c2 = sH < sL;
                c3 = abs(mL) <= cfg.criteria.nearMargin_psi && abs(mH) <= cfg.criteria.nearMargin_psi;
                cr(end+1, :) = {fr.name, g, T0.BIP_CO2_C1(1), T0.BIP_CO2_nC8(1), models(m).name, ...
                    T0.xCO2_exactMatch(o(1)), T0.xCO2_exactMatch(o(end)), sL, sH, mL, mH, c1, c2, c3, c1 && c2 && c3}; %#ok<AGROW>
            end
        end
    end
    frontCrit = textToString(cell2table(cr, 'VariableNames', {'Front', 'Point', 'BIP_CO2_C1', 'BIP_CO2_nC8', ...
        'Model', 'xCO2_low_exact', 'xCO2_high_exact', 'Shift_low_psi', 'Shift_high_psi', 'MissLow_psi', ...
        'MissHigh_psi', 'BothPositive', 'HighLessThanLow', 'NearBands', 'AllPass'}));

    % BIP-only end of each front: nominal composition, BIP solved per mixture (bracketed on the grid)
    bo = {};
    for fr = fronts
        for qa = 1:numel(iCO2pos)
            a = iCO2pos(qa);  c = dp.getCase(caseIdx(a));
            gN = PnomGrid.(fr.name)(qa, :) - ref(a).bulkExp;
            j = find(isfinite(gN(1:end-1)) & isfinite(gN(2:end)) & gN(1:end-1) .* gN(2:end) <= 0, 1);
            val = NaN;  b1 = NaN;  b8 = NaN;
            if ~isempty(j)
                fun = @(v) bulkWithBIP(cfg, fr.name, v, c, dp, caseIdx) - ref(a).bulkExp;
                try
                    val = fzero(fun, fr.grid([j, j+1]), optimset('TolX', 1e-4));
                catch
                    val = fr.grid(j) - gN(j) * (fr.grid(j+1) - fr.grid(j)) / (gN(j+1) - gN(j));
                end
                [b1, b8] = frontBIPs(fr.name, val, cfg);
            end
            bo(end+1, :) = {fr.name, c.Mixture, c.z(find(c.Components == "CO2", 1)), b1, b8, ~isempty(j)}; %#ok<AGROW>
        end
    end
    bipOnly = textToString(cell2table(bo, 'VariableNames', {'Front', 'Mixture', 'xCO2_nominal', ...
        'BIP_CO2_C1_required', 'BIP_CO2_nC8_required', 'BracketedOnGrid'}));

    fprintf('\n=== Stage E: exact-match fronts (bulk matched exactly at every point) ===\n');
    disp(frontPts);  disp(frontCrit);  disp(bipOnly);
end

%% 8. Write ---------------------------------------------------------------------------------
modelTable = struct2table(rmfield(models, {'rp', 'xRange', 'rpManual'}));
modelTable.rp_mid_nm = arrayfun(@(s) s.rp.mid, models).';
modelTable.rp_lo_nm  = arrayfun(@(s) s.rp.lo,  models).';
modelTable.rp_hi_nm  = arrayfun(@(s) s.rp.hi,  models).';
out = fullfile(cfg.resultsDir, cfg.outFile);
save(out + ".mat", 'cfg', 'models', 'runTable', 'critTable', 'tipTable', 'tcTable', 'bulkModel', 'xEff', ...
    'frontPts', 'frontRuns', 'frontCrit', 'bipOnly');
xls = out + ".xlsx";
if isfile(xls), delete(xls); end
writetable(modelTable, xls, 'Sheet', 'Models');
if ~isempty(tcTable), writetable(tcTable, xls, 'Sheet', 'A_CO2_Tc_blind'); end
writetable(runTable,  xls, 'Sheet', 'BC_runs');
writetable(critTable, xls, 'Sheet', 'D_criteria');
writetable(tipTable,  xls, 'Sheet', 'Pair_TIPs_used');
if ~isempty(frontPts),  writetable(frontPts,  xls, 'Sheet', 'E_front_points'); end
if ~isempty(frontRuns), writetable(frontRuns, xls, 'Sheet', 'E_front_runs'); end
if ~isempty(frontCrit), writetable(frontCrit, xls, 'Sheet', 'E_front_criteria'); end
if ~isempty(bipOnly),   writetable(bipOnly,   xls, 'Sheet', 'E_BIP_only_end'); end
fprintf('Written: %s\n', xls);

%% ============================== Local functions ==============================

function rp = loadRp(mm, cfg)
% r_p [nm] for the model. cfg.rpSource = "manual": rpManual (shift basis). "file": the
% back-calculation file first, manual values as fallback.
rp = struct('mid', NaN, 'lo', NaN, 'hi', NaN, 'src', "none");
if cfg.rpSource == "manual"
    if isempty(mm.rpManual), return; end
    rp = struct('mid', mm.rpManual(1), 'lo', mm.rpManual(2), 'hi', mm.rpManual(3), 'src', "manual (shift basis)");
    return;
end
fn = fullfile(cfg.resultsDir, mm.rpFile);
if isfile(fn)
    S = load(fn, 'derived');
    key = matlab.lang.makeValidName(mm.form);
    if isfield(S, 'derived') && isfield(S.derived, key)
        d = S.derived.(key);
        rp = struct('mid', 1e9 * d.rp_mid_m, 'lo', 1e9 * d.rp_lo_m, 'hi', 1e9 * d.rp_hi_m, 'src', "file " + mm.rpFile);
    end
end
if ~isempty(mm.rpManual)
    if isfinite(rp.mid) && abs(rp.mid - mm.rpManual(1)) > 0.02
        warning('co2:RpMismatch', '%s: file r_p %.3f nm differs from manual %.3f nm; the file is used.', mm.name, rp.mid, mm.rpManual(1));
    elseif ~isfinite(rp.mid)
        rp = struct('mid', mm.rpManual(1), 'lo', mm.rpManual(2), 'hi', mm.rpManual(3), 'src', "manual fallback");
    end
end
end

function co = tipCoefficients(tab, comb, mode, form, rock)
sel = string(tab.Combination) == comb & string(tab.Mode) == mode & string(tab.Form) == form;
if nnz(sel) ~= 1
    error('co2:CorrelationRow', '%d rows match %s / %s / %s in the Correlation sheet.', nnz(sel), comb, mode, form);
end
vn = "alpha_" + rock;
if ~ismember(vn, string(tab.Properties.VariableNames))
    error('co2:NoRockOffset', 'Column %s not found: no offset for rock %s.', vn, rock);
end
co = struct('alpha', tab.(vn)(sel), 'B', tab.B(sel), 'descriptor', string(tab.Descriptor(sel)));
end

function X = pairDescriptor(fluid, name, T, kB)
w  = fluid.omega(:);
e  = fluid.LJ_Energy(:);
MW = fluid.MW(:);
Ts = (kB * T) ./ sqrt(e * e.');
switch name
    case "wbar",      X = 0.5 * (w + w.');
    case "inv_wbar",  X = 2 ./ (w + w.');
    case "Tstar",     X = Ts;
    case "inv_Tstar", X = 1 ./ Ts;
    case "ln_Tstar",  X = log(Ts);
    case "inv_Mgeo",  X = 1 ./ sqrt(MW * MW.');
    otherwise
        error('co2:UnknownDescriptor', 'Descriptor "%s" is not implemented.', name);
end
end

function [K, T] = buildTIP(fluid, co, Tk, kB, xRange, clip)
% TIP_ij = -exp(alpha_rock + B * x_ij); clip = true limits x to the calibrated range.
X  = pairDescriptor(fluid, co.descriptor, Tk, kB);
Xu = X;
if clip, Xu = min(max(X, xRange(1)), xRange(2)); end
n  = size(X, 1);
K  = -exp(co.alpha + co.B * Xu);
K(1:n+1:end) = 0;
K  = (K + K.') / 2;
nm = string(fluid.ComponentNames);
rows = cell(0, 5);
for a = 1:n-1
    for b = a+1:n
        flag = "within";
        if X(a, b) < xRange(1) - 1e-9, flag = "below"; end
        if X(a, b) > xRange(2) + 1e-9, flag = "above"; end
        rows(end+1, :) = {nm(a) + "-" + nm(b), X(a, b), Xu(a, b), flag, K(a, b)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Pair', 'x', 'x_used', 'RangeFlag', 'TIP'});
T.Pair = string(T.Pair);  T.RangeFlag = string(T.RangeFlag);
end

function [inBand, miss] = bandMiss(s, lo, up)
% Signed distance outside the band [lo, up] (0 inside; negative below, positive above).
inBand = s >= lo && s <= up;
if s < lo, miss = s - lo; elseif s > up, miss = s - up; else, miss = 0; end
end

function x = extractXco2(zStr)
% CO2 mole fraction from the "z" strings ("0.2000,0.7200,0.0800"; CO2 is the first component).
x = zeros(numel(zStr), 1);
for q = 1:numel(zStr)
    v = str2double(split(zStr(q), ","));
    x(q) = v(1);
end
end

function [xr, scanX, scanP] = effectiveCO2(runner, c, Pexp_psi, cfg)
% Every CO2 fraction (C1:nC8 ratio held) whose model bulk Pdew equals the experimental bulk.
psi = entities.DewPointDataset.psi2Pa;
iC = find(c.Components == "CO2", 1);
rest = c.z;  rest(iC) = 0;  rest = rest / sum(rest);
scanX = cfg.composition.scanGrid(:).';
scanP = nan(size(scanX));
guess = c.PdewBulk;
for k = 1:numel(scanX)
    [scanP(k), ok] = bulkAtX(runner, c, rest, iC, scanX(k), guess, cfg.composition.xRound);
    if ok, guess = scanP(k) * psi; end
end
g = scanP - Pexp_psi;
xr = [];
for j = 1:numel(scanX) - 1
    if isfinite(g(j)) && isfinite(g(j+1)) && g(j) * g(j+1) < 0
        a = scanX(j);  b = scanX(j+1);
        fun = @(x) bulkAtX(runner, c, rest, iC, x, c.PdewBulk, cfg.composition.xRound) - Pexp_psi;
        try
            x = fzero(fun, [a b], optimset('TolX', 2e-4));
        catch
            x = a - g(j) * (b - a) / (g(j+1) - g(j));
        end
        xr(end+1) = round(x, cfg.composition.xRound); %#ok<AGROW>
    end
end
end

function [P_psi, ok] = bulkAtX(runner, c, rest, iC, x, guess_Pa, nd)
psi = entities.DewPointDataset.psi2Pa;
x = round(x, nd);
z = (1 - x) * rest;  z(iC) = x;
c2 = c;  c2.z = z(:);  c2.PdewBulk = guess_Pa;
[Pb, ~, ok] = runner.bulk(c2);
P_psi = NaN;
if ok && isfinite(Pb), P_psi = Pb / psi; else, ok = false; end
end

function [b1, b8] = frontBIPs(front, v, cfg)
% BIP pair (CO2-C1, CO2-nC8) at front coordinate v.
switch front
    case "CO2C1", b1 = v;                              b8 = cfg.front.bipLit.CO2_nC8;
    case "both",  b1 = cfg.front.bipLit.CO2_C1 + v;    b8 = cfg.front.bipLit.CO2_nC8 + v;
    otherwise, error('co2:Front', 'Unknown front %s.', front);
end
end

function writeBipWorkbook(src, dst, changes)
% Copy the mixture workbook and set symmetric BIPs in sheet "BIP" (row/column labels = components).
copyfile(src, dst, 'f');
C = readcell(src, 'Sheet', 'BIP');
miss = cellfun(@(e) isa(e, 'missing'), C);
C(miss) = {''};
hdr = string(C(1, 2:end));  lab = string(C(2:end, 1));
for q = 1:size(changes, 1)
    a = string(changes{q, 1});  b = string(changes{q, 2});  v = changes{q, 3};
    ra = find(lab == a, 1);  rb = find(lab == b, 1);  ca = find(hdr == a, 1);  cb = find(hdr == b, 1);
    if any(cellfun(@isempty, {ra, rb, ca, cb}))
        error('co2:BipLabel', 'BIP sheet has no row/column for %s or %s.', a, b);
    end
    C{1 + ra, 1 + cb} = v;  C{1 + rb, 1 + ca} = v;
end
writecell(C, dst, 'Sheet', 'BIP', 'Range', 'A1');
end

function P_psi = bulkWithBIP(cfg, front, v, c, dp, caseIdx)
% Model bulk Pdew of case c at nominal composition with the front BIPs at coordinate v.
psi = entities.DewPointDataset.psi2Pa;
[b1, b8] = frontBIPs(front, v, cfg);
f = fullfile(cfg.front.workDir, "bip_bracket.xlsx");
writeBipWorkbook(cfg.mixtureFile, f, {"CO2", "C1", b1; "CO2", "nC8", b8});
rr = solvers.DewPointRunner(f, cfg.poreRadiusNominal, 'maxShift_psi', Inf);
rr.engine(c.Mixture, c.Rock);
[Pb, ~, ok] = rr.bulk(c);
P_psi = NaN;  if ok, P_psi = Pb / psi; end
end

function t = textToString(t)
for v = string(t.Properties.VariableNames)
    col = t.(v);
    if iscell(col) && ~isempty(col) && all(cellfun(@(e) isstring(e) || ischar(e), col))
        t.(v) = string(col);
    end
end
end