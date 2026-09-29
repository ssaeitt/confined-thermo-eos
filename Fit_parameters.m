%% FIT_PARAMETERS  Global calibration of the confined EOS (OOP port of V31.0)
%
% Model: a_eff = a_bulk - c  |  a_mix_eff = a_mix_bulk - c_mix  (thermo.ConfinedEOS)
%
% Decision vector (V31 ordering)
%   x = [k, p, lambda, TIP_FWIOnly(1..G), TIP_Combined(1..G)]
%   G binary (pair, rock) groups, default:
%   C1-nC5|EF2, C1-nC8|EF2, C1-nC10|EF2, C1-nC8|B1
%
% Objective (V31 form):
%   J = RMSE(dTc) / wTc + mean_over_modes( RMSE(Pdew) [psi] ) / wPd
% k, p, lambda are shared by both modes (mode-independent); each mode has its
% own TIP set (TIP = the kijc property of ConfinedEOS).
%
% cfg.strategy = "joint"      : all parameters by fmincon/SQP (V31)
%              = "sequential" : k, p, lambda on dTc only, then 1-D TIP fits
% cfg.tcModel  = "linearized" : V31 expression dTc = c_i(Tc,conf,exp) / a_i(Tc)
%              = "exact"      : confined PR critical point (ConfinedEOS.pureCriticalPoint)
% Both dTc measures are always reported for the final parameters.

clear; clc; close all;

%% Preflight: locate +entities/+thermo/+solvers and put their parent on the path
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end          % run by selection / section
if isfolder(fullfile(scriptDir, "src")), addpath(fullfile(scriptDir, "src")); end
for p = ["entities", "thermo", "solvers"]
    if isempty(meta.package.fromName(p))
        cand = [dir(fullfile(scriptDir, "+" + p)); dir(fullfile(scriptDir, "**", "+" + p))];
        cand = unique(string({cand.folder}));
        cand = cand(endsWith(cand, filesep + "+" + p));
        if isempty(cand)
            error('fit:PackageNotFound', ['Folder "+%s" not found under %s.\n' ...
                'Move Fit_parameters.m to the folder that contains +entities, +thermo, +solvers.'], ...
                p, scriptDir);
        end
        addpath(fileparts(cand(1)));
        fprintf('Added to path: %s  (for +%s)\n', fileparts(cand(1)), p);
    end
end
required = ["entities.FluidProperties", "entities.RockProperties", "entities.TcShiftDataset", ...
            "entities.DewPointDataset", "thermo.ConfinedEOS", "solvers.StabilityTester", ...
            "solvers.FlashEngine"];
missing = required(arrayfun(@(c) isempty(meta.class.fromName(c)), required));
if ~isempty(missing)
    error('fit:MissingClasses', ['Package folders found, but these classes are missing or ' ...
        'fail to parse: %s.\nCheck the file names match the class names exactly.'], ...
        strjoin(missing, ', '));
end

%% 0. Configuration ---------------------------------------------------------
cfg = struct();
% Project layout:  <root>/Fit_parameters.m
%                  <root>/src/+entities, +thermo, +solvers
%                  <root>/config/MixtureData.xlsx, Extended_Dataset.xlsx
%                  <root>/results/            (created if absent)
cfg.configDir    = fullfile(scriptDir, "config");
cfg.resultsDir   = fullfile(scriptDir, "results");
cfg.mixtureFile  = fullfile(cfg.configDir, "MixtureData.xlsx");       % replaces "Data Single.xlsx"
cfg.tcShiftFile  = fullfile(cfg.configDir, "Extended_Dataset - Original.xlsx");
cfg.tcShiftSheet = "All Data";

% CALIBRATION WALL ONLY. This pseudo-mineral wall is used for the literature
% dTc data (silica pores) and for nothing else: every mixture engine keeps the
% real per-mineral mineralogy of its own rock (EF2 / B1), loaded by
% RockProperties, and the confined dew points use that.
cfg.wallSourceRock = "EF2";                     % row used only to read mineral_E
cfg.wallMinerals   = "Quartz";                  % quartz-only calibration wall
                                                % (V31 used ["Quartz","Plagioclase","K-feldspar"])

cfg.T          = 293.15;
cfg.poreRadius = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});   % [m], reported radii

cfg.fitPairs = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"];
cfg.fwiStart = [64.44, 1.3745, 10.0338];                  % [k, p, lambda]
cfg.fwiLB    = [0.1, 0.1,  0.0];
cfg.fwiUB    = [300, 2.5, 20.0];
cfg.tipStart = [-3.1759, -2.4635, -1.7946, -0.7284];      % same start for both modes

% Stage 1 (dTc) numerics
cfg.tcLambdaSeeds   = [0 2 4.5 7 10];        % multi-start seeds in lambda (at k = 100 and 150)
cfg.tcLambdaProfile = [0 1 2 3 4 5 6 8 10];  % RMSE(lambda) landscape with k, p re-optimised; [] to skip
cfg.tcFDStep        = 1e-6;                  % explicit finite-difference step for Stage-1 gradients
cfg.tipLB    = -100;
cfg.tipUB    = 50;

cfg.modes    = ["FWIOnly", "Combined"];
cfg.strategy = "joint";
cfg.tcModel  = "linearized";
cfg.wTc      = 0.05;
cfg.wPd      = 100;                   % [psi]
cfg.failurePenalty_psi = 1000;        % NaN = skip failed points (V31 behaviour)
cfg.pdTarget     = "absolute";        % "absolute": Pmodel - Pexp,mid (V31)
                                      % "shift"   : (Pmodel - Pbulk,model) - (Pexp,mid - Pbulk,exp)
cfg.combinedFromFWI = true;           % seed Combined solves from the FWI-only solution
cfg.ssFallback      = true;           % on a Newton failure: SS -> Newton -> Broyden -> cold Newton
cfg.maxShift_psi    = 500;            % |Pconf - Pbulk| above this is a lost branch, not a solution
cfg.contSteps       = [0.25 0.5 0.75 1];   % ConfinementScale homotopy, started from the bulk state
cfg.sumXTol         = 1e-6;           % |sum(z/K) - 1| at the returned root
cfg.diagGroups      = ["C1-nC5|EF2"];  % groups whose per-case c_mix diagnostics are printed
cfg.pGuessFactor = 1.1;               % cold-start guess = factor * bulk Pdew (V31)
cfg.flashTol     = 1e-10;             % tight inner solves for finite-difference gradients
cfg.fdStep       = 1e-6;              % fmincon relative FD step
cfg.maxIter      = 200;
cfg.tipScanPoints = 31;               % sequential strategy only
cfg.makePlots    = true;
cfg.runTag       = "";                % set below from strategy + tcModel
cfg.outFile      = "fit_results";     % "_<runTag>" is appended automatically
cfg.runLog       = "fit_runs_log.csv";% one appended row per run, for cross-run comparison


% --- Run identity: outputs are tagged, so a stale file can never be mistaken
%     for a new result and the three combinations cannot overwrite each other.
if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
cfg.runTag  = cfg.strategy + "-" + cfg.tcModel;
cfg.outFile = fullfile(cfg.resultsDir, cfg.outFile + "_" + cfg.runTag);
scriptsSeen = which('Fit_parameters', '-all');
fprintf('\n================= RUN: %s =================\n', cfg.runTag);
fprintf('  strategy = %s | tcModel = %s\n', cfg.strategy, cfg.tcModel);
fprintf('  script   : %s\n', string(scriptsSeen{1}));
if numel(scriptsSeen) > 1
    warning('fit:DuplicateScript', ['Fit_parameters exists in %d places on the path; the one above ' ...
        'is the one that ran:\n  %s'], numel(scriptsSeen), strjoin(string(scriptsSeen), "\n  "));
end
fprintf('  data     : %s\n             %s\n', cfg.mixtureFile, cfg.tcShiftFile);
fprintf('  output   : %s.[mat|xlsx]\n', cfg.outFile);
if isfile(cfg.outFile + ".mat")
    d = dir(cfg.outFile + ".mat");
    fprintf(2, '  NOTE: overwriting the existing result for this tag (written %s)\n', d.date);
end

psi = entities.DewPointDataset.psi2Pa;

%% 1. Entities ----------------------------------------------------------------
fluidAll = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile);
rockRef  = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallSourceRock);

% Calibration wall: pseudo-mineral with the arithmetic-mean silicate energy.
% With a single mineral, ConfinedEOS gives eps_fw = sqrt(eps_f * E_wall), as in V31.
Em = arrayfun(@(m) getfield(rockRef.mineral(m), 'Epsilon_K'), cfg.wallMinerals);     %#ok<GFLD>
Rm = arrayfun(@(m) getfield(rockRef.mineral(m), 'GrainDensity'), cfg.wallMinerals);  %#ok<GFLD>
calRock = entities.RockProperties("CalWall", 0, "CalWall", mean(Em), mean(Rm), 1);
eosCal  = thermo.ConfinedEOS(fluidAll, calRock);   % Stage 1 only
fprintf('Calibration wall (dTc only): %s, E = %.4f K. Mixture engines keep their own mineralogy.\n', ...
        strjoin(cfg.wallMinerals, "+"), mean(Em));

tcAll = entities.TcShiftDataset.loadFromWorkbook(cfg.tcShiftFile, 'Sheet', cfg.tcShiftSheet);
tc    = tcAll.select(ismember(tcAll.Component, fluidAll.ComponentNames));
[tcIdx, tcR] = tc.modelInputs(fluidAll);          % r = sigma_model / d_H (d_H used as tabulated)

dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
groupKey = dp.Mixture + "|" + dp.Rock;
isFit = dp.HasConfined & dp.NumComponents == 2 & ismember(groupKey, cfg.fitPairs);
isVal = dp.HasConfined & dp.NumComponents > 2;
unlisted = setdiff(unique(groupKey(dp.HasConfined & dp.NumComponents == 2)), cfg.fitPairs);
if ~isempty(unlisted)
    warning('fit:UnlistedPairs', 'Binary data not in cfg.fitPairs (ignored): %s', strjoin(unlisted, ", "));
end
for g = cfg.fitPairs
    if ~any(isFit & groupKey == g)
        error('fit:NoData', 'No confined binary data for %s.', g);
    end
end

G = numel(cfg.fitPairs);
M = numel(cfg.modes);
fitCases = find(isFit);
caseGroup = arrayfun(@(i) find(cfg.fitPairs == groupKey(i)), fitCases);

% Engines: one ConfinedEOS / StabilityTester / FlashEngine per (mixture, rock)
engines = containers.Map();
for i = [fitCases; find(isVal)].'
    key = char(groupKey(i));
    if ~isKey(engines, key)
        fl = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile, dp.Mixture(i));
        rk = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, dp.Rock(i));
        eos = thermo.ConfinedEOS(fl, rk);
        st  = solvers.StabilityTester(eos);
        fe  = solvers.FlashEngine(eos, st, 'Tolerance', cfg.flashTol);
        engines(key) = struct('eos', eos, 'flash', fe);
    end
end

% Bulk dew points (r = Inf: c = 0, independent of k, p, lambda and TIP)
bulkModelPa = nan(dp.NumCases, 1);
bulkModelK  = cell(dp.NumCases, 1);
for i = [fitCases; find(isVal)].'
    c = dp.getCase(i);
    eng = engines(char(groupKey(i)));
    eng.eos.kijc = zeros(numel(c.z));
    eng.eos.ConfinementScale = 1;
    [Pb, Kb, ~, ~, st] = eng.flash.solveDewPoint(c.T, c.PdewBulk, c.z, Inf, 'CapMode', "FWIOnly");
    if st.converged
        bulkModelPa(i) = Pb;
        bulkModelK{i}  = Kb;
    else
        warning('fit:BulkFailed', 'Bulk dew point not converged for %s (%s): %s', ...
                dp.Mixture(i), dp.Rock(i), st.reason);
    end
end
if cfg.pdTarget == "shift" && any(isnan(bulkModelPa(fitCases)))
    error('fit:BulkFailed', 'Bulk dew point failed for a fitted case; "shift" target unavailable.');
end

ctx = struct('cfg', cfg, 'bulkModelPa', bulkModelPa, 'bulkModelK', {bulkModelK}, 'eosCal', eosCal, 'tc', tc, 'tcIdx', tcIdx, 'tcR', tcR, ...
             'dp', dp, 'fitCases', fitCases, 'caseGroup', caseGroup, 'groupKey', groupKey, ...
             'engines', engines, 'warm', containers.Map(), 'G', G, 'M', M);

% Provenance: inputs that Stage 1 depends on, so runs can be compared exactly
dM = dir(cfg.mixtureFile);  dT = dir(cfg.tcShiftFile);
provenance = struct('when', string(datetime('now')), 'strategy', cfg.strategy, ...
    'tcModel', cfg.tcModel, 'mixtureFile', cfg.mixtureFile, 'tcShiftFile', cfg.tcShiftFile, ...
    'mixtureBytes', dM.bytes, 'mixtureDate', string(dM.date), ...
    'tcShiftBytes', dT.bytes, 'tcShiftDate', string(dT.date), ...
    'nTcPoints', tc.NumPoints, 'sum_dH', sum(tc.dH_model), 'sum_dTc', sum(tc.dTc), ...
    'calWallEnergy_K', mean(Em), 'calMinerals', strjoin(cfg.wallMinerals, "+"), ...
    'LJ_Size', fluidAll.LJ_Size, 'LJ_Energy', fluidAll.LJ_Energy, ...
    'fwiStart', cfg.fwiStart, 'fwiLB', cfg.fwiLB, 'fwiUB', cfg.fwiUB, ...
    'tcFDStep', cfg.tcFDStep, 'matlab', string(version));
fprintf(['\nProvenance: %s (%d bytes, %s) | %s (%d bytes, %s)\n' ...
         '            %d dTc points, sum(d_H) = %.4f, sum(dTc) = %.4f, wall E = %.4f K\n'], ...
        cfg.mixtureFile, dM.bytes, dM.date, cfg.tcShiftFile, dT.bytes, dT.date, ...
        tc.NumPoints, sum(tc.dH_model), sum(tc.dTc), mean(Em));

%% 2. Optimisation --------------------------------------------------------------
x0 = [cfg.fwiStart, repmat(cfg.tipStart, 1, M)];
lb = [cfg.fwiLB,    repmat(cfg.tipLB, 1, G * M)];
ub = [cfg.fwiUB,    repmat(cfg.tipUB, 1, G * M)];

tcProfTable = table();
scanEdge = false(G, M);           % sequential only: scan minimum at a bound
scanFlat = false(G, M);           % sequential only: penalty-flat objective
switch cfg.strategy
    case "joint"
        fprintf('\nJoint optimisation: %d parameters, %d dTc points, %d Pdew points x %d modes\n', ...
            numel(x0), tc.NumPoints, numel(fitCases), M);
        opts = optimoptions('fmincon', 'Display', 'iter-detailed', 'Algorithm', 'sqp', ...
            'MaxIterations', cfg.maxIter, 'FiniteDifferenceStepSize', cfg.fdStep);
        [xBest, fBest, exitflag] = fmincon(@(x) jointObjective(x, ctx), x0, [], [], [], [], ...
                                           lb, ub, [], opts);

    case "sequential"
        fprintf('\nStage 1: k, p, lambda on dTc (%s)\n', cfg.tcModel);
        fwi = fitTcOnly(ctx);
        if ~isempty(cfg.tcLambdaProfile)
            tcProfTable = tcLambdaProfile(ctx);
            fprintf('  RMSE(lambda) with k, p re-optimised (Stage-1 landscape):\n');
            disp(tcProfTable);
            [~, iB] = min(tcProfTable.RMSE);
            if tcProfTable.RMSE(iB) < tcRMSE(fwi, ctx) - 1e-6
                fprintf(2, ['  Stage 1 stopped above the profile minimum ' ...
                            '(RMSE %.5f at lambda = %.2f); adopting the profile point.\n'], ...
                        tcProfTable.RMSE(iB), tcProfTable.lambda(iB));
                fwi = [tcProfTable.k(iB), tcProfTable.p(iB), tcProfTable.lambda(iB)];
            end
            setFWI(ctx, fwi);
        end
        setFWI(ctx, fwi);
        tips = zeros(G, M);
        for m = 1:M
            for g = 1:G
                fobj = @(t) groupSSE(t, g, cfg.modes(m), ctx);
                grid = linspace(cfg.tipLB, cfg.tipUB, cfg.tipScanPoints);
                vals = arrayfun(fobj, grid);
                [~, j] = min(vals);
                scanEdge(g, m) = (j == 1 || j == numel(grid));
                scanFlat(g, m) = max(vals) - min(vals) <= 1e-9 * max(1, abs(min(vals)));
                tips(g, m) = fminbnd(fobj, grid(max(j-1, 1)), grid(min(j+1, end)), optimset('TolX', 1e-6));
                fprintf('  %-12s %-9s TIP = %+.5f%s\n', cfg.fitPairs(g), cfg.modes(m), tips(g, m), ...
                        string(repmat(' (objective flat: no accepted points)', 1, scanFlat(g, m))));
            end
        end
        xBest = [fwi, tips(:).'];
        fBest = jointObjective(xBest, ctx);
        exitflag = NaN;

    otherwise
        error('fit:Strategy', 'Unknown strategy "%s".', cfg.strategy);
end

%% 3. Final evaluation --------------------------------------------------------
[fwi, tips] = unpack(xBest, ctx);
[~, parts] = jointObjective(xBest, ctx);

dTcLin = predictTc(ctx, "linearized");
dTcEx  = predictTc(ctx, "exact");
tcTable = [tc.toTable(), table(dTcLin, dTcEx, 'VariableNames', {'dTc_linearized', 'dTc_exactEOS'})];
tcStats = @(p) [sqrt(mean((p - tc.dTc).^2, 'omitnan')), 100 * mean(abs(p - tc.dTc) ./ tc.dTc, 'omitnan')];
sLin = tcStats(dTcLin);  sEx = tcStats(dTcEx);

caseRows = {};
for m = 1:M
    [res, det] = pdResiduals(ctx, cfg.modes(m), tips(:, m));  %#ok<ASGLU>
    for n = 1:numel(det)
        d = det(n);
        c = dp.getCase(d.case);
        dg = struct('cm_V', NaN, 'cm_L', NaN, 'ratio_V', NaN, 'ratio_L', NaN);
        if d.ok
            dg = diagnoseCase(engines(char(groupKey(d.case))), c, cfg.poreRadius(char(c.Rock)), ...
                              d.Pdew, exp(d.lnK), d.Pcap, cfg.modes(m));
        end
        caseRows(end+1, :) = {cfg.modes(m), dp.Mixture(d.case), dp.Rock(d.case), ...
            "[" + strjoin(string(dp.MoleFrac{d.case}.'), ",") + "]", ...
            dp.MoleFrac{d.case}(end), tips(caseGroup(n), m), ...
            dp.BulkPdew_psi(d.case), dp.ConfMid_psi(d.case), dp.ConfHalfWidth_psi(d.case), d.Pdew / psi, ...
            dp.ShiftMid_psi(d.case), d.Pdew / psi - dp.BulkPdew_psi(d.case), d.Pcap / psi, ...
            abs(d.Pdew / psi - dp.ConfMid_psi(d.case)) <= dp.ConfHalfWidth_psi(d.case), d.ok, ...
            d.route, d.reason, d.sumX - 1, d.Z_V, d.Z_L, d.rhoM_V, d.rhoM_L, ...
            dg.cm_V, dg.cm_L, dg.ratio_V, dg.ratio_L}; %#ok<SAGROW>
    end
end
caseTable = cell2table(caseRows, 'VariableNames', {'Mode', 'Mixture', 'Rock', 'z', 'y_heavy', 'TIP', ...
    'Bulk_exp_psi', 'Pdew_exp_psi', 'HalfWidth_psi', 'Pdew_model_psi', 'Shift_exp_psi', 'Shift_model_psi', ...
    'Pcap_psi', 'InBand', 'Accepted', 'Route', 'Reason', 'sumX_res', 'Z_V', 'Z_L', ...
    'rhoM_V', 'rhoM_L', 'cm_V', 'cm_L', 'cm_over_am_V', 'cm_over_am_L'});
caseTable = convertvars(caseTable, {'Mode', 'Mixture', 'Rock', 'z', 'Route', 'Reason'}, 'string');

tipTable = table(cfg.fitPairs.', accumarray(caseGroup(:), 1, [G, 1]), 'VariableNames', {'Group', 'nPoints'});
tipsReport = tips;
for m = 1:M
    md = cfg.modes(m);
    sel = caseTable.Mode == md;
    err = caseTable.Pdew_model_psi(sel) - caseTable.Pdew_exp_psi(sel);
    rmseG = sqrt(accumarray(caseGroup(:), err.^2, [G, 1], @(v) mean(v, 'omitnan')));
    nAcc = accumarray(caseGroup(:), double(caseTable.Accepted(sel)), [G, 1]);
    nPts = accumarray(caseGroup(:), 1, [G, 1]);
    tolB = 1e-3 * (cfg.tipUB - cfg.tipLB);
    atB  = abs(tips(:, m) - cfg.tipLB) <= tolB | abs(tips(:, m) - cfg.tipUB) <= tolB | scanEdge(:, m);
    status = repmat("ok", G, 1);
    status(nAcc < nPts) = "partial (" + string(nAcc(nAcc < nPts)) + "/" + string(nPts(nAcc < nPts)) + " accepted)";
    status(atB & nAcc > 0) = status(atB & nAcc > 0) + ", at bound";
    noFit = nAcc == 0 | scanFlat(:, m);
    status(noFit) = "no fit (0 accepted points)";
    tipOut = tips(:, m);
    tipOut(noFit) = NaN;
    tipsReport(:, m) = tipOut;
    rmseG(noFit) = NaN;
    tipTable.("TIP_" + md)       = tipOut;
    tipTable.("Status_" + md)    = status;
    tipTable.("nAccepted_" + md) = nAcc;
    tipTable.("RMSE_psi_" + md)  = rmseG;
    tipTable.("AtBound_" + md)   = atB & ~noFit;
end

% Bulk model check
bulkModel = bulkModelPa / psi;
bulkTable = table(dp.Mixture, dp.Rock, dp.BulkPdew_psi, bulkModel, bulkModel - dp.BulkPdew_psi, ...
    'VariableNames', {'Mixture', 'Rock', 'Bulk_exp_psi', 'Bulk_model_psi', 'Dev_psi'});

% Multicomponent validation with fitted binary TIPs
valRows = {};
for i = find(isVal).'
    c = dp.getCase(i);
    eng = engines(char(groupKey(i)));
    for m = 1:M
        [K, note] = assembleTIP(cfg.fitPairs, tipsReport(:, m), c.Components, c.Rock);
        P = NaN; Pc = NaN;
        if ~isempty(K)
            eng.eos.kijc = K;
            [P, ~, Pc, ~, st] = eng.flash.solveDewPoint(c.T, cfg.pGuessFactor * c.PdewBulk, c.z, ...
                                    cfg.poreRadius(char(c.Rock)), 'CapMode', cfg.modes(m));
            if ~st.converged, P = NaN; end
        end
        valRows(end+1, :) = {cfg.modes(m), c.Mixture, c.Rock, dp.ConfMid_psi(i), P / psi, ...
            P / psi - dp.ConfMid_psi(i), Pc / psi, note}; %#ok<SAGROW>
    end
end
valTable = cell2table(valRows, 'VariableNames', {'Mode', 'Mixture', 'Rock', 'Pdew_exp_psi', ...
    'Pdew_model_psi', 'Dev_psi', 'Pcap_psi', 'Note'});
valTable = convertvars(valTable, {'Mode', 'Mixture', 'Rock', 'Note'}, 'string');

%% 4. Report ----------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf(' GLOBAL OPTIMAL PARAMETERS  (strategy: %s, dTc model: %s, exitflag %g)\n', ...
        cfg.strategy, cfg.tcModel, exitflag);
fprintf('  k = %.4f | p = %.4f | lambda = %.4f\n', fwi);
fprintf('  J = %.5f  [RMSE_dTc = %.4f; RMSE_Pdew = %s psi]\n', fBest, parts.rmseTc, ...
        strjoin(compose("%.2f (%s)", parts.rmsePd(:), cfg.modes(:)), ", "));
fprintf('  Failed dew points: %s\n', strjoin(compose("%d/%d (%s)", ...
        arrayfun(@(md) sum(~caseTable.Accepted(caseTable.Mode == md)), cfg.modes(:)), ...
        numel(fitCases) * ones(M, 1), cfg.modes(:)), ", "));
fprintf('  dTc with final parameters: linearized RMSE %.4f / AARD %.1f %%;  exact-EOS RMSE %.4f / AARD %.1f %%\n', ...
        sLin, sEx);
fprintf('------------------------------------------------------------------\n');
disp(tipTable);
fprintf('Per-case results:\n');   disp(caseTable);
for gd = cfg.diagGroups
    pr = split(gd, "|");
    sel = caseTable(caseTable.Mixture == pr(1) & caseTable.Rock == pr(2), :);
    if isempty(sel), continue; end
    fprintf('\nc_mix diagnostics for %s (cm = confinement term of the mixture, am = bulk attraction):\n', gd);
    disp(sel(:, {'Mode', 'y_heavy', 'TIP', 'Pdew_model_psi', 'Shift_model_psi', 'Accepted', 'Route', ...
                 'cm_V', 'cm_L', 'cm_over_am_V', 'cm_over_am_L', 'Z_V', 'Z_L', 'rhoM_V', 'rhoM_L', ...
                 'sumX_res', 'Reason'}));
end
fprintf('Bulk check:\n');         disp(bulkTable);
fprintf('Multicomponent validation:\n'); disp(valTable);

fwiTable = table(["k"; "p"; "lambda"], fwi(:), 'VariableNames', {'Parameter', 'Value'});
save(cfg.outFile + ".mat", 'cfg', 'provenance', 'tcProfTable', 'xBest', 'fBest', 'fwi', 'tips', 'tipsReport', 'parts', 'tcTable', ...
     'tipTable', 'caseTable', 'bulkTable', 'valTable');
xls = cfg.outFile + ".xlsx";
writetable(fwiTable,  xls, 'Sheet', 'FWI_params');
writetable(tipTable,  xls, 'Sheet', 'TIP');
writetable(caseTable, xls, 'Sheet', 'Pdew_cases');
writetable(tcTable,   xls, 'Sheet', 'dTc');
writetable(bulkTable, xls, 'Sheet', 'Bulk_check');
writetable(valTable,  xls, 'Sheet', 'Validation');
writetable(struct2table(structfun(@(v) string(mat2str(v, 6)), provenance, 'UniformOutput', false)), ...
           xls, 'Sheet', 'Provenance');
if height(tcProfTable) > 0, writetable(tcProfTable, xls, 'Sheet', 'dTc_lambda_profile'); end

if cfg.makePlots
    try
        plotTc(tc, dTcLin, dTcEx, fwi, cfg);
        plotPdew(caseTable, cfg);
    catch ME
        warning('fit:PlotFailed', 'Plotting failed (%s); continuing with Stage 6.', ME.message);
    end
end
fprintf('\nMain results written: %s.mat / .xlsx\n', cfg.outFile);

% --- cross-run log: one row per run, so two runs that produce identical
%     parameters under different tags are impossible to miss
logFile = fullfile(cfg.resultsDir, cfg.runLog);
logRow = table(string(datetime('now')), cfg.runTag, cfg.strategy, cfg.tcModel, ...
    fwi(1), fwi(2), fwi(3), 'VariableNames', {'When', 'RunTag', 'Strategy', 'tcModel', ...
    'k', 'p', 'lambda'});
for mm = 1:M
    for gg = 1:G
        logRow.(matlab.lang.makeValidName("TIP_" + cfg.modes(mm) + "_" + cfg.fitPairs(gg))) = tips(gg, mm);
    end
end
if isfile(logFile)
    prev = readtable(logFile, 'TextType', 'string');
    same = abs(prev.k - fwi(1)) < 1e-8 & abs(prev.p - fwi(2)) < 1e-8 & abs(prev.lambda - fwi(3)) < 1e-8;
    if any(same & prev.RunTag ~= cfg.runTag)
        warning('fit:IdenticalAcrossRuns', ['This run produced the same k, p, lambda as run tag(s) ' ...
            '%s. Different strategy/tcModel settings must not give identical parameters - check that ' ...
            'cfg.strategy and cfg.tcModel were actually changed before this run.'], ...
            strjoin(unique(prev.RunTag(same & prev.RunTag ~= cfg.runTag)), ", "));
    end
    logRow = [prev; logRow];
end
writetable(logRow, logFile);
fprintf('Run log: %s (%d runs recorded)\n', logFile, height(logRow));
fprintf('  %s: k = %.6f, p = %.6f, lambda = %.6f\n', cfg.runTag, fwi);
fprintf('  TIPs: %s\n', strjoin(compose("%s %s = %+.5f", ...
        repmat(cfg.fitPairs(:), M, 1), repelem(cfg.modes(:), G, 1), tips(:)), " | "));

%% ============================== Local functions ==============================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end














function [fwi, tips] = unpack(x, ctx)
    fwi  = x(1:3);                                   % [k, p, lambda]
    tips = reshape(x(4:end), ctx.G, ctx.M);          % columns follow cfg.modes
end

function setFWI(ctx, fwi)
    % V31 order [k, p, lambda] -> ConfinedEOS properties
    eosList = [{ctx.eosCal}, cellfun(@(e) e.eos, values(ctx.engines), 'UniformOutput', false)];
    for n = 1:numel(eosList)
        eosList{n}.k       = fwi(1);
        eosList{n}.pT_wall = fwi(2);
        eosList{n}.lambda  = fwi(3);
    end
end

function [J, parts] = jointObjective(x, ctx)
    [fwi, tips] = unpack(x, ctx);
    setFWI(ctx, fwi);
    rTc = predictTc(ctx, ctx.cfg.tcModel) - ctx.tc.dTc;
    rTc(~isfinite(rTc)) = 1;                          % failed critical point
    rmseTc = sqrt(mean(rTc.^2));
    rmsePd = zeros(1, ctx.M);
    for m = 1:ctx.M
        [r, d] = pdResiduals(ctx, ctx.cfg.modes(m), tips(:, m));
        rmsePd(m) = sqrt(objValue(r, d, ctx.cfg) / numel(r));
    end
    J = rmseTc / ctx.cfg.wTc + mean(rmsePd) / ctx.cfg.wPd;
    parts = struct('rmseTc', rmseTc, 'rmsePd', rmsePd);
end

function fwi = fitTcOnly(ctx)
    % Least-squares fit of [k, p, lambda] to dTc. Multi-start lsqnonlin;
    % the result is checked never to be worse than the start point.
    cfg = ctx.cfg;
    resFun = @(f) tcResiduals(f, ctx);
    f0 = tcRMSE(cfg.fwiStart, ctx);
    lamSeeds = cfg.tcLambdaSeeds(:);
    starts = [cfg.fwiStart; [100 + 0 * lamSeeds, 1.3 + 0 * lamSeeds, lamSeeds]; ...
              [150 + 0 * lamSeeds, 1.3 + 0 * lamSeeds, lamSeeds]];
    best = f0;  fwi = cfg.fwiStart;
    for s = 1:size(starts, 1)
        if exist('lsqnonlin', 'file') == 2
            % Central differences with an explicit step: the default step
            % (~1.5e-8 relative) is below the fzero resolution of the confined
            % critical point, which makes the gradient noisy and the result
            % run-dependent.
            o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, ...
                             'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
                             'MaxFunctionEvaluations', 5000, 'FiniteDifferenceType', 'central', ...
                             'FiniteDifferenceStepSize', cfg.tcFDStep);
            x = lsqnonlin(resFun, starts(s, :), cfg.fwiLB, cfg.fwiUB, o);
        else
            o = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
            x = fmincon(@(f) tcRMSE(f, ctx), starts(s, :), [], [], [], [], cfg.fwiLB, cfg.fwiUB, [], o);
        end
        fx = tcRMSE(x, ctx);
        fprintf('  start %d: RMSE %.5f -> %.5f   [k p lambda] = [%.4f %.4f %.4f]\n', ...
                s, tcRMSE(starts(s, :), ctx), fx, x);
        if fx < best, best = fx;  fwi = x; end
    end
    r2 = tcR2(fwi, ctx);
    fprintf('  Stage 1: RMSE %.5f -> %.5f, R2 = %.4f, [k p lambda] = [%.4f %.4f %.4f]\n', ...
            f0, best, r2, fwi);
end

function r2 = tcR2(fwi, ctx)
    setFWI(ctx, fwi);
    p = predictTc(ctx, ctx.cfg.tcModel);
    ok = isfinite(p);
    r2 = 1 - sum((p(ok) - ctx.tc.dTc(ok)).^2) / sum((ctx.tc.dTc(ok) - mean(ctx.tc.dTc(ok))).^2);
end

function T = tcLambdaProfile(ctx)
    % RMSE(lambda) with k and p re-optimised at each lambda. Shows whether a
    % Stage-1 result sits at the global minimum or is a premature stop.
    cfg = ctx.cfg;
    lam = cfg.tcLambdaProfile(:);
    n = numel(lam);
    [kk, pp, rr, R2] = deal(nan(n, 1));
    for i = 1:n
        f = @(v) tcResiduals([v(1), v(2), lam(i)], ctx);
        o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, ...
                         'StepTolerance', 1e-12, 'FiniteDifferenceType', 'central', ...
                         'FiniteDifferenceStepSize', cfg.tcFDStep);
        v = lsqnonlin(f, [150, 1.3], cfg.fwiLB(1:2), cfg.fwiUB(1:2), o);
        kk(i) = v(1);  pp(i) = v(2);
        rr(i) = tcRMSE([v(1), v(2), lam(i)], ctx);
        R2(i) = tcR2([v(1), v(2), lam(i)], ctx);
    end
    T = table(lam, kk, pp, rr, R2, 'VariableNames', {'lambda', 'k', 'p', 'RMSE', 'R2'});
end

function r = tcResiduals(fwi, ctx)
    setFWI(ctx, fwi);
    r = predictTc(ctx, ctx.cfg.tcModel) - ctx.tc.dTc;
    r(~isfinite(r)) = 1;
end

function f = tcRMSE(fwi, ctx)
    setFWI(ctx, fwi);
    r = predictTc(ctx, ctx.cfg.tcModel) - ctx.tc.dTc;
    r(~isfinite(r)) = 1;
    f = sqrt(mean(r.^2));
end

function pred = predictTc(ctx, model)
    eos = ctx.eosCal;
    pred = nan(ctx.tc.NumPoints, 1);
    for n = 1:ctx.tc.NumPoints
        i  = ctx.tcIdx(n);
        Tc = eos.Fluid.Tc(i);
        try
            switch model
                case "linearized"   % V31: c evaluated at the experimental confined Tc
                    [ac, ~, ~] = eos.pureParameters(Tc, Inf);
                    [~, ~, c]  = eos.pureParameters(Tc * (1 - ctx.tc.dTc(n)), ctx.tcR(n));
                    pred(n) = c(i) / ac(i);
                case "exact"
                    pred(n) = 1 - eos.pureCriticalPoint(i, ctx.tcR(n)) / Tc;
            end
        catch ME
            if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
        end
    end
end

function f = groupSSE(t, g, mode, ctx)
    tips = zeros(ctx.G, 1);
    tips(g) = t;
    [r, d] = pdResiduals(ctx, mode, tips, find(ctx.caseGroup == g));
    f = objValue(r, d, ctx.cfg);
end

function f = objValue(res, det, cfg)
    % Optimiser objective only: accepted SSE plus a fixed penalty per rejected
    % point, so the search is pushed away from regions with lost branches.
    % Statistics (RSS, RMSE, F-test) must use accepted points only - never this.
    f = sum(res.^2, 'omitnan') + cfg.failurePenalty_psi^2 * sum(~[det.ok]);
end

function st = fitStats(res, det)
    % Accepted-only statistics for one (mode, beta) fit.
    ok = [det.ok].';
    st = struct('rss', sum(res(ok).^2, 'omitnan'), 'N', nnz(ok), ...
                'nRejected', nnz(~ok), 'acceptedCases', [det(ok).case], ...
                'rmse', sqrt(mean(res(ok).^2, 'omitnan')));
end

function [res, det] = pdResiduals(ctx, mode, tips, subset)
    % Residual [psi] = Pdew_model - Pdew_exp,mid. V31 formed
    % (Pmodel - Pbulk,exp) - (Pmid - Pbulk,exp), which is identical.
    if nargin < 4, subset = 1:numel(ctx.fitCases); end
    psi = entities.DewPointDataset.psi2Pa;
    res = nan(numel(subset), 1);
    det = struct('case', {}, 'Pdew', {}, 'Pcap', {}, 'ok', {}, 'route', {}, 'reason', {}, ...
                 'sumX', {}, 'Z_V', {}, 'Z_L', {}, 'rhoM_V', {}, 'rhoM_L', {}, 'lnK', {}, ...
                 'cm_V', {}, 'cm_L', {}, 'group', {}, 'y_heavy', {});
    for n = 1:numel(subset)
        i   = ctx.fitCases(subset(n));
        g   = ctx.caseGroup(subset(n));
        c   = ctx.dp.getCase(i);
        eng = ctx.engines(char(ctx.groupKey(i)));
        eng.eos.kijc = [0, tips(g); tips(g), 0];
        r   = ctx.cfg.poreRadius(char(c.Rock));
        bulk = struct('P', ctx.bulkModelPa(i), 'K', ctx.bulkModelK{i});
        if ~isfinite(bulk.P), bulk.P = c.PdewBulk; bulk.K = []; end
        [P, Pcap, ok, sInfo] = solveConfined(eng, c, r, mode, ctx.warm, ...
                                             sprintf('%d|%s', i, mode), ctx.cfg, bulk);
        cmV = NaN;  cmL = NaN;
        if ok      % c_mix contrast: depends on composition, T and r only
            x   = c.z ./ exp(sInfo.lnK);
            cmV = eng.eos.confinementMix(c.z, c.T, r);
            cmL = eng.eos.confinementMix(x / sum(x), c.T, r);
        end
        if ok
            if ctx.cfg.pdTarget == "shift"
                res(n) = ((P - ctx.bulkModelPa(i)) - (c.PdewConfMid - c.PdewBulk)) / psi;
            else
                res(n) = (P - c.PdewConfMid) / psi;
            end
        else
            res(n) = NaN;        % rejected: never a pseudo-residual (see objValue)
        end
        det(n) = struct('case', i, 'Pdew', P, 'Pcap', Pcap, 'ok', ok, ...
                        'route', sInfo.route, 'reason', string(sInfo.reason), ...
                        'sumX', sInfo.sumX, 'Z_V', sInfo.Z_V, 'Z_L', sInfo.Z_L, ...
                        'rhoM_V', sInfo.rhoM_V, 'rhoM_L', sInfo.rhoM_L, 'lnK', sInfo.lnK, ...
                        'cm_V', cmV, 'cm_L', cmL, 'group', g, 'y_heavy', c.z(end));
    end
end

function [P, Pcap, ok, info] = solveConfined(eng, c, r, mode, warm, key, cfg, bulk)
    % Returns an ACCEPTED dew point only (acceptRoot). Routes, first success wins:
    %   1) warm      : this case's last accepted state (same mode);
    %   2) fwi-cont. : Combined only, seeded by the FWI-only solution at the same TIP;
    %   3) direct    : full c, seeded by the converged BULK state (K, P);
    %   4) homotopy  : fallback only, ConfinementScale stepped up from the bulk
    %                  state (the bulk problem is not re-solved or re-validated).
    ok = false;  P = NaN;  Pcap = NaN;
    info = emptyInfo("none");
    args = {'CapMode', mode, 'Solver', "newton"};
    eng.eos.ConfinementScale = 1;
    reasons = strings(0, 1);

    if isKey(warm, key)
        w = warm(key);
        [P, K, Pcap, st] = tryFlash(eng, c, r, w.P, w.K, w.Pcap, args, cfg);
        [ok, info] = acceptRoot(P, K, st, c, bulk.P, cfg);
        info.route = "warm";
        if ~ok, reasons(end+1) = "warm: " + info.reason; end
    end

    if ~ok && mode == "Combined" && cfg.combinedFromFWI
        fwiKey = char(replace(string(key), "Combined", "FWIOnly"));
        if isKey(warm, fwiKey)
            w = warm(fwiKey);
            [P, K, Pcap, st] = tryFlash(eng, c, r, w.P, w.K, [], args, cfg);
            [ok, info] = acceptRoot(P, K, st, c, bulk.P, cfg);
            info.route = "fwi-continuation";
            if ~ok, reasons(end+1) = "fwi-continuation: " + info.reason; end
        end
    end

    if ~ok
        [P, K, Pcap, st] = tryFlash(eng, c, r, cfg.pGuessFactor * bulk.P, bulk.K, [], args, cfg);
        [ok, info] = acceptRoot(P, K, st, c, bulk.P, cfg);
        info.route = "direct";
        if ~ok, reasons(end+1) = "direct: " + info.reason; end
    end

    if ~ok
        [P, Pcap, ok, info] = homotopy(eng, c, r, cfg, bulk, args);
        if ~ok, reasons(end+1) = info.reason; end
    end

    eng.eos.ConfinementScale = 1;
    if ok
        warm(key) = struct('P', P, 'K', exp(info.lnK), 'Pcap', Pcap);
    else
        P = NaN;  Pcap = NaN;
        info.reason = strjoin(reasons, " | ");
    end
end

function info = emptyInfo(route)
    info = struct('route', route, 'reason', "", 'sumX', NaN, 'normLnK', NaN, ...
                  'Z_V', NaN, 'Z_L', NaN, 'rhoM_V', NaN, 'rhoM_L', NaN, 'shift_psi', NaN, 'lnK', []);
end

function [P, Pcap, ok, info] = homotopy(eng, c, r, cfg, bulk, args)
    % Continuation in ConfinementScale starting FROM the converged bulk state.
    % s = 0 is solved only if no bulk state exists.
    P = NaN;  Pcap = NaN;  ok = false;
    info = emptyInfo("homotopy");
    steps = cfg.contSteps(cfg.contSteps > 0);
    if isempty(bulk.K)
        steps = [0, steps];
        Kprev = [];
    else
        Kprev = bulk.K;
    end
    Pprev = bulk.P;  Pcprev = [];
    for s = steps
        eng.eos.ConfinementScale = s;
        [Ps, Ks, Pcs, st] = tryFlash(eng, c, r, Pprev, Kprev, Pcprev, args, cfg);
        [okS, infoS] = acceptRoot(Ps, Ks, st, c, bulk.P, cfg);
        if ~okS
            eng.eos.ConfinementScale = 1;
            info.reason = "homotopy: failed at ConfinementScale = " + string(s) + " (" + infoS.reason + ")";
            return;
        end
        Kprev = Ks;  Pprev = Ps;  Pcprev = Pcs;
        [P, Pcap, info] = deal(Ps, Pcs, infoS);
        info.route = "homotopy";
    end
    eng.eos.ConfinementScale = 1;
    ok = true;
end

function [P, K, Pcap, st] = tryFlash(eng, c, r, Pguess, Kseed, Pcapseed, args, cfg)
    % Newton from the given seeds; on failure, fall back in this order:
    %   (a) successive substitution from the same seeds, then Newton from the SS
    %       state (SS has no Jacobian and no line search, so it survives the
    %       regions where Newton stalls with "line-search-failure");
    %   (b) Broyden from the SS state (cheaper Jacobian refreshes);
    %   (c) Newton with the engine's own stability-test seeding.
    % Fallbacks only supply better seeds: acceptRoot still judges the result.
    if nargin < 8, cfg = struct('ssFallback', true); end
    extra = {};
    if ~isempty(Kseed),    extra = [extra, {'K_seed', Kseed, 'PreconditionSS', false}]; end
    if ~isempty(Pcapseed), extra = [extra, {'Pcap_seed', Pcapseed}]; end
    [P, K, Pcap, ~, st] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, args{:}, extra{:});
    if st.converged || ~isfield(cfg, 'ssFallback') || ~cfg.ssFallback
        return;
    end

    % (a) successive substitution as a seed generator (its own tolerance is
    % tight, so it usually reports not-converged; the state is still useful)
    ssArgs = replaceSolver(args, "ss");
    [Pss, Kss, Pcss, ~, stss] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, ssArgs{:}, extra{:});
    if ~isfinite(Pss) || any(~isfinite(Kss)) || norm(log(max(Kss, 1e-14))) < 0.01
        Pss = Pguess;  Kss = Kseed;  Pcss = Pcapseed;      % SS unusable: keep the original seeds
    end
    ex2 = {'K_seed', Kss, 'PreconditionSS', false};
    if ~isempty(Pcss) && isfinite(Pcss), ex2 = [ex2, {'Pcap_seed', Pcss}]; end
    if ~isempty(Kss)
        [P2, K2, Pc2, ~, st2] = eng.flash.solveDewPoint(c.T, Pss, c.z, r, args{:}, ex2{:});
        if st2.converged || betterState(st2, st)
            [P, K, Pcap, st] = deal(P2, K2, Pc2, st2);
        end
        if st.converged, return; end

        % (b) Broyden from the same SS state
        qnArgs = replaceSolver(args, "quasinewton");
        [P3, K3, Pc3, ~, st3] = eng.flash.solveDewPoint(c.T, Pss, c.z, r, qnArgs{:}, ex2{:});
        if st3.converged || betterState(st3, st)
            [P, K, Pcap, st] = deal(P3, K3, Pc3, st3);
        end
        if st.converged, return; end
    end

    % (c) cold Newton with stability-test seeding
    [P4, K4, Pc4, ~, st4] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, args{:});
    if st4.converged || betterState(st4, st)
        [P, K, Pcap, st] = deal(P4, K4, Pc4, st4);
    end
end

function a = replaceSolver(args, solver)
    a = args;
    i = find(strcmp(a, 'Solver'), 1);
    if isempty(i)
        a = [a, {'Solver', solver}];
    else
        a{i + 1} = solver;
    end
end

function tf = betterState(stNew, stOld)
    % Prefer the state with the smaller residual norm when neither converged.
    tf = isfield(stNew, 'residual_norm') && isfield(stOld, 'residual_norm') && ...
         isfinite(stNew.residual_norm) && ~(stOld.residual_norm <= stNew.residual_norm);
end

function [ok, info] = acceptRoot(P, K, st, c, Pbulk, cfg)
    % A converged residual is not sufficient: the solver can land on another
    % root. Require dew closure, phase identity (FlashEngine.PhaseIdentityTest:
    % mass density / heavy enrichment, never Z) and a physically possible shift.
    psi = entities.DewPointDataset.psi2Pa;
    info = emptyInfo("");
    if isempty(K) || any(~isfinite(K))
        info.reason = "solver: " + string(st.reason);
        ok = false;  return;
    end
    sumX = sum(c.z ./ K);
    info.sumX = sumX;  info.normLnK = norm(log(K));  info.lnK = log(K);
    info.Z_V = st.Z_V;  info.Z_L = st.Z_L;
    if isfield(st, 'rhoM_V'), info.rhoM_V = st.rhoM_V;  info.rhoM_L = st.rhoM_L; end
    info.shift_psi = (P - Pbulk) / psi;
    ok = false;
    if ~st.converged
        info.reason = "solver: " + string(st.reason);
    elseif ~isfinite(P) || P <= 0
        info.reason = "non-finite P";
    elseif abs(sumX - 1) > cfg.sumXTol
        info.reason = sprintf("sum(z/K) - 1 = %.2e", sumX - 1);
    elseif ~st.incipientIsLiquid
        info.reason = sprintf("phase identity (rho_m,L = %.1f vs rho_m,V = %.1f kg/m3)", ...
                              info.rhoM_L, info.rhoM_V);
    elseif abs(info.shift_psi) > cfg.maxShift_psi
        info.reason = sprintf("|shift| = %.0f psi > %.0f psi (lost branch)", ...
                              abs(info.shift_psi), cfg.maxShift_psi);
    else
        ok = true;
        info.reason = "accepted";
    end
end

function d = diagnoseCase(eng, c, r, P, K, Pcap, mode)
    % c_mix of both phases at the reported root (amix_c from ConfinedEOS).
    x = c.z ./ K;  x = x / sum(x);
    PL = P;
    if mode == "Combined", PL = max(P - Pcap, 1e3); end
    [~, ~, ~, pv] = eng.eos.calculateState(P,  c.T, c.z, -1, r);
    [~, ~, ~, pl] = eng.eos.calculateState(PL, c.T, x,   +1, r);
    d = struct('cm_V', pv.amix_c, 'cm_L', pl.amix_c, 'amix_V', pv.amix_a, 'amix_L', pl.amix_a, ...
               'ratio_V', pv.amix_c / pv.amix_a, 'ratio_L', pl.amix_c / pl.amix_a);
end

function [K, note] = assembleTIP(fitPairs, tips, comps, rock)
    % kijc for a multicomponent case from fitted binary TIPs of the same rock.
    % Pairs with C1 must have been fitted; other pairs default to 0.
    n = numel(comps);
    K = zeros(n);
    note = "";
    for a = 1:n-1
        for b = a+1:n
            g = find(fitPairs == comps(a) + "-" + comps(b) + "|" + rock | ...
                     fitPairs == comps(b) + "-" + comps(a) + "|" + rock, 1);
            if ~isempty(g) && isnan(tips(g))
                K = [];
                note = "TIP for " + comps(a) + "-" + comps(b) + " in " + rock + " is 'no fit'";
                return;
            elseif ~isempty(g)
                K(a, b) = tips(g);  K(b, a) = tips(g);
            elseif any([comps(a), comps(b)] == "C1")
                K = [];
                note = "no fitted TIP for " + comps(a) + "-" + comps(b) + " in " + rock;
                return;
            else
                note = strtrim(note + " " + comps(a) + "-" + comps(b) + "=0 (assumed)");
            end
        end
    end
end

function plotTc(tc, dTcLin, dTcEx, fwi, cfg)
    figure('Name', 'Confinement: dTc', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 19 13]);
    hold on; grid on; box on;
    comps = unique(tc.Component, 'stable');
    col = lines(numel(comps));
    for i = 1:numel(comps)
        s = tc.Component == comps(i);
        plot(tc.dH(s), tc.dTc(s), 'o', 'Color', col(i, :), 'MarkerSize', 4, 'DisplayName', comps(i) + " Exp");
        plot(tc.dH(s), dTcLin(s), 'x', 'Color', col(i, :), 'MarkerSize', 9, 'LineWidth', 1.5, ...
             'DisplayName', comps(i) + " linearized");
        plot(tc.dH(s), dTcEx(s), '+', 'Color', col(i, :), 'MarkerSize', 7, 'DisplayName', comps(i) + " exact EOS");
    end
    legend('Location', 'northwest', 'NumColumns', 3);
    xlabel('\sigma_{LJ}/r_p');  ylabel('\Delta T_c');
    title(sprintf('Pure-component suppression (%s fit): k=%.2f, p=%.3f, \\lambda=%.2f', cfg.tcModel, fwi));
end

function plotPdew(caseTable, cfg)
    figure('Name', 'Confinement: mixture Pdew', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    styles = ["r-o", "b-s", "g-^"];
    for g = cfg.fitPairs
        parts = split(g, "|");
        nexttile; hold on; grid on; box on;
        yAll = [];
        for m = 1:numel(cfg.modes)
            sel = caseTable(caseTable.Mixture == parts(1) & caseTable.Rock == parts(2) & ...
                            caseTable.Mode == cfg.modes(m), :);
            if isempty(sel), continue; end
            [yH, o] = sort(sel.y_heavy);
            sel = sel(o, :);
            if m == 1
                errorbar(yH, sel.Shift_exp_psi, sel.HalfWidth_psi, 'ks', ...
                         'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Exp.');
                yAll = [yAll; sel.Shift_exp_psi + sel.HalfWidth_psi; sel.Shift_exp_psi - sel.HalfWidth_psi];
            end
            ok = isfinite(sel.Shift_model_psi);
            plot(yH(ok), sel.Shift_model_psi(ok), styles(m), 'LineWidth', 1.5, 'DisplayName', cfg.modes(m));
            if any(~ok)   % mark non-converged cases on the axis
                plot(yH(~ok), zeros(nnz(~ok), 1), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5, ...
                     'DisplayName', cfg.modes(m) + " (failed)");
            end
            yAll = [yAll; sel.Shift_model_psi(ok)];
        end
        title(parts(1) + " (" + parts(2) + ")");
        xlabel('y_{heavy}');  ylabel('\Delta P_{dew} (psi)');
        yAll = yAll(isfinite(yAll));
        if ~isempty(yAll)
            pad = max(0.08 * range(yAll), 1);
            ylim([min(yAll) - pad, max(yAll) + pad]);
        end
        xl = xlim;  xlim(xl + 0.02 * diff(xl) * [-1 1]);
        if g == cfg.fitPairs(1), legend('Location', 'best'); end
    end
    sgtitle('Mixture dew-point shift: FWI-only vs Combined', 'FontWeight', 'bold');
end

function full = resolveFile(name, root)
    % Absolute path of a data file: as given if it exists, otherwise the unique
    % match under root (recursive).
    name = string(name);
    if isfile(name)
        info = dir(name);
        full = string(fullfile(info.folder, info.name));
        return;
    end
    [~, base, ext] = fileparts(name);
    hits = dir(fullfile(root, "**", base + ext));
    hits = hits(~[hits.isdir]);
    if isempty(hits)
        error('fit:FileNotFound', '"%s" not found under %s.', base + ext, root);
    elseif numel(hits) > 1
        error('fit:AmbiguousFile', '"%s" found in several folders:\n  %s\nSet cfg with the full path.', ...
            base + ext, strjoin(string({hits.folder}), "\n  "));
    end
    full = string(fullfile(hits.folder, hits.name));
end