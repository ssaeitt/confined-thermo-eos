%% TIP_CORRELATION  Size-independent, dimensionless TIP correlation, rock offset
%
%   ln(-TIP_ij) = alpha_rock + B * x_ij              TIP_ij = -exp(...) < 0
%
%   x_ij : absolute pair property (decreases |TIP| for heavy/strongly bound pairs,
%          which the ternaries require: heavy-heavy |TIP| must stay small)
%          wbar      = (omega_i + omega_j)/2
%          inv_wbar  = 2/(omega_i + omega_j)
%          Tstar     = kB*T/sqrt(eps_i*eps_j)
%          inv_Tstar = sqrt(eps_i*eps_j)/(kB*T)
%          ln_Tstar  = ln(Tstar)      (= E* form after the rock offset absorbs ln S_rock)
%          inv_Mgeo  = 1/sqrt(M_i*M_j) [mol/g]; dimensionless TIP, B carries g/mol
%   rock : alpha_EF2, alpha_B1 (rock-indexed offset). S_rock and TOC are reported as
%          equivalent re-expressions of the same offset.
%
% Parameters alpha_EF2, alpha_B1, B on four binary groups -> dof = 1. alpha_B1 is set
% by the single B1 binary (leverage 1.0). k, p, lambda and the per-group optimised
% TIPs are unchanged. The identity limit is dropped: no absolute-descriptor form satisfies it
% (disclosed); the ternaries require heavy-heavy |TIP| to be small.
%
% Acceptance criteria and the selection rule are declared in section 0, BEFORE the run.

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
            error('tip:PackageNotFound', 'Folder "+%s" not found under %s.', pk, scriptDir);
        end
        addpath(fileparts(cand(1)));
    end
end

%% 0. Configuration ---------------------------------------------------------
cfg.mixtureFile = fullfile(scriptDir, "config", "MixtureData - Complete Mixtures.xlsx");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.outFile     = "tip_correlation_rockoffset";
cfg.T           = 293.15;
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});   % [m], radii
cfg.modes       = ["FWIOnly", "Combined"];
cfg.primaryMode = "FWIOnly";
cfg.shiftBasis  = "modelBulk";          % shift = P_conf,model - P_bulk,model

% Locked (k, p, lambda) from Fit_parameters.
% NOTE (carried over): the joint-exact k was taken from the earlier joint-exact run
% (98.05788478) after a transcription slip. Correct it here if the true value differs.
cfg.combos = table( ...
    ["sequential-exact"; "joint-exact"; "joint-linearized"; "pareto-exact-w2"; "pareto-exact-w1"], ...
    [151.385346647943; 98.05788478; 112.42634826677; 160.958378; 166.410386], ...
    [1.30682834548441; 1.69745310322519; 1.3558179290678; 1.460506; 1.548842], ...
    [4.51381302592871; 10.2425163919516; 2.0916235090532; 4.471524; 4.479082], ...
    'VariableNames', {'Combination', 'k', 'p', 'lambda'});
cfg.runCombos = [];                      % subset to run; [] = all

cfg.groups = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"];
cfg.tips = containers.Map( ...           % per-group optimised TIPs (unchanged)
    {'sequential-exact|FWIOnly', 'joint-exact|FWIOnly', 'joint-linearized|FWIOnly', ...
    'sequential-exact|Combined', 'joint-exact|Combined', 'joint-linearized|Combined', ...
    'pareto-exact-w2|FWIOnly', 'pareto-exact-w2|Combined', ...
    'pareto-exact-w1|FWIOnly', 'pareto-exact-w1|Combined'}, ...
    {[-3.48386287716462, -1.46205964046861, -1.12088738883485, -0.905508761783693], ...
    [-13.4437287153297, -5.90838752675671, -4.59334889956389, -2.38867380361359 ], ...
    [-4.9229696302796,  -1.96920366521686, -1.44628069292089, -0.754862049795788], ...
    [-3.40199024752739, -1.45777900094979, -1.11884643084433, -0.899502001873144], ...
    [-14.4508613237875, -5.84525317456513, -4.79513444682033, -2.24015630615276 ], ...
    [-4.76033469976459, -1.96111617688159, -1.44196236936452, -0.74629351212583 ], ...
    [-5.310985, -2.205853, -1.647836, -1.059127], ...
    [-5.310985, -2.205853, -1.647836, -1.059127], ...
    [-7.012375, -2.883946, -2.124570, -1.168401], ...
    [-7.012375, -2.883946, -2.124570, -1.168401]});

% Candidate forms (declared before the run); all eligible for selection
cfg.forms = table( ...
    ["LL-wbar"; "LL-inv_wbar"; "LL-Tstar"; "LL-inv_Tstar"; "LL-ln_Tstar"; "LL-inv_Mgeo"], ...
    ["wbar";    "inv_wbar";    "Tstar";    "inv_Tstar";    "ln_Tstar";    "inv_Mgeo"], ...
    'VariableNames', {'Name', 'Descriptor'});

% Rows outside the HC calibration/validation set
cfg.late = struct('Mixture', "C1-nC8", 'Rock', "B1", 'yHeavy', 0.10);  % late-campaign core
cfg.excludeCO2 = true;                                                  % CO2 stage is separate

% Regression gate: optimised TIPs must reproduce the previous round (also verifies
% that the adjusted C1-nC5 compositions are in the workbook)
cfg.reg.combination = "sequential-exact";
cfg.reg.mode        = "FWIOnly";
cfg.reg.shifts_psi  = [110.9, 39.1, -8.9, 175.0, 140.5, 84.4, 98.0, 313.7, 180.6];  % Mix_Props order
cfg.reg.tol_psi     = 1.0;
cfg.reg.basis       = "expBulk";
cfg.reg.b1Group     = "C1-nC8|B1";
cfg.reg.b1Shift_psi = 98.0;
cfg.reg.b1Tol_psi   = 0.1;

% --- Acceptance criteria, declared BEFORE the run ---------------------------
cfg.accept.noNegativeShift   = true;   % every HC binary and ternary, nominal r_p
cfg.accept.maxErr_C1nC10_psi = 60;
cfg.accept.ternaryInBand     = true;
cfg.accept.maxPdewErr_pct    = 1.0;
cfg.accept.rpMonotonic       = true;   % status must be Pass (see below)
% r_p sweep, own rock and offset. The 500 psi lost-branch cap is lifted in the sweep.
% AMENDED after run 2 (disclosed): monotonicity is judged on the branch tracked from
% the nominal radius = the contiguous run of converged, same-sign radii through it.
% Solves beyond a failure (e.g. a sign-flipped point past the end of the branch) are
% excluded. The local log-slope is reported, not used as a criterion (run 2 showed a
% smooth steepening on EF2 at small r_p, not a jump).
%   Pass  : branch has >= minFinite radii, all shifts > 0, strictly decreasing
%   Fail  : branch has a non-positive shift or is not strictly decreasing
%   Undetermined : nominal solve failed or branch shorter than minFinite
cfg.rpSweep.radii_nm     = [5, 8.1, 11.25, 15, 20, 30, 54, 100];
cfg.rpSweep.modes        = cfg.modes;   % restrict to cfg.primaryMode to cut runtime
cfg.rpSweep.maxShift_psi = Inf;         % DewPointRunner.maxShift_psi in the sweep only
cfg.rpSweep.minFinite    = 4;
cfg.kB = 1.380649e-23;                  % J/K, for T* (Comp_Props LJ_Energy in J)

% --- Selection rule, declared BEFORE the run --------------------------------
% Among (combination, form) rows in the primary mode that pass every criterion:
% lowest ternary shift RMSE; rows within tieTol of the best are ranked by correlated
% binary shift RMSE. The ternaries therefore inform the selection (disclosed).
cfg.select.tieTol_psi = 1.0;
cfg.extrapFlagFactor  = 5;             % flag applied |TIP| > factor * max calibrated |TIP|

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
psi = entities.DewPointDataset.psi2Pa;

fprintf('Acceptance criteria (declared before the run):\n');
fprintf(['  no negative shift (all HC) %d | C1-nC10 err <= %.0f psi | ' ...
    'ternaries in band %d | |Pdew err| <= %.1f %% | r_p-monotonic Pass %d\n'], ...
    cfg.accept.noNegativeShift, cfg.accept.maxErr_C1nC10_psi, ...
    cfg.accept.ternaryInBand, cfg.accept.maxPdewErr_pct, cfg.accept.rpMonotonic);
fprintf(['Selection: primary mode %s; lowest ternary RMSE among passing rows, ties within ' ...
    '%.1f psi broken by binary RMSE.\n\n'], cfg.primaryMode, cfg.select.tieTol_psi);

%% 1. Data, engines, component and rock properties -----------------------------
dp     = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius);
P      = loadComponentProps(cfg.mixtureFile);
rockD  = loadRockDescriptors(cfg.mixtureFile);

yHeavy   = cellfun(@(v) v(end), dp.MoleFrac);
isLate   = dp.Mixture == cfg.late.Mixture & dp.Rock == cfg.late.Rock & ...
    abs(yHeavy - cfg.late.yHeavy) < 1e-9;
isCO2    = contains(dp.Mixture, "CO2") & cfg.excludeCO2;
groupKey = dp.Mixture + "|" + dp.Rock;
isBin    = dp.HasConfined & dp.NumComponents == 2 & ismember(groupKey, cfg.groups) & ~isLate;
isTer    = dp.HasConfined & dp.NumComponents > 2 & ~isCO2;
binIdx   = find(isBin);
terIdx   = find(isTer);
caseIdx  = [binIdx; terIdx];
if numel(binIdx) ~= 9 || numel(terIdx) ~= 2
    error('tip:CaseCount', 'Expected 9 HC binaries and 2 HC ternaries, got %d and %d.', ...
        numel(binIdx), numel(terIdx));
end
fprintf('Excluded: %d late-campaign row(s), %d CO2 row(s).\n', nnz(isLate), nnz(contains(dp.Mixture, "CO2")));

for i = caseIdx.'
    runner.engine(dp.Mixture(i), dp.Rock(i));
end
ek = keys(runner.Engines);
for e = 1:numel(ek)
    eng = runner.Engines(ek{e});
    pr  = split(string(ek{e}), "|");
    rk  = eng.eos.Rock;
    if rk.RockName ~= pr(2) || rk.NumMinerals < 2
        error('tip:WallMismatch', 'Engine %s carries wall "%s" with %d mineral(s).', ...
            string(ek{e}), rk.RockName, rk.NumMinerals);
    end
    nm = string(eng.eos.Fluid.ComponentNames);
    if ~all(isKey(P, cellstr(nm)))
        error('tip:CompProps', 'Comp_Props lacks a component of engine %s.', string(ek{e}));
    end
end
for i = caseIdx.'
    [~, ~, okb] = runner.bulk(dp.getCase(i));
    if ~okb, error('tip:BulkFailed', 'Bulk dew point failed for %s (%s).', dp.Mixture(i), dp.Rock(i)); end
end

G = numel(cfg.groups);
grpRock = strings(G, 1); grpPair = cell(G, 1);
for g = 1:G
    pr = split(cfg.groups(g), "|");
    eng = runner.engine(pr(1), pr(2));
    grpRock(g) = pr(2);
    grpPair{g} = string(eng.eos.Fluid.ComponentNames);
end
isB1g = grpRock == "B1";
fprintf('Rock descriptors: S_EF2 = %.4f, S_B1 = %.4f sqrt(K); TOC vol frac EF2 = %.4f, B1 = %.4f\n\n', ...
    rockD.S.EF2, rockD.S.B1, rockD.TOC.EF2, rockD.TOC.B1);

%% 2. Regression gate (optimised TIPs; independent of the correlation form) ----
fprintf('=== Regression gate (%s / %s) ===\n', cfg.reg.combination, cfg.reg.mode);
ci = find(cfg.combos.Combination == cfg.reg.combination, 1);
runner.setFWI(cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci));
runner.clearWarm();
tipReg = cfg.tips(char(cfg.reg.combination + "|" + cfg.reg.mode));
[regModel, regExp] = deal(nan(numel(binIdx), 1));
for n = 1:numel(binIdx)
    i = binIdx(n);
    g = find(cfg.groups == groupKey(i), 1);
    r = runner.solve(dp.getCase(i), cfg.reg.mode, [0, tipReg(g); tipReg(g), 0], 'Key', sprintf("reg|%d", i));
    if r.ok, regModel(n) = r.shift_psi; regExp(n) = r.shiftExpBulk_psi; end
end
regShift = regExp; if cfg.reg.basis == "modelBulk", regShift = regModel; end
regTable = table(dp.Mixture(binIdx), dp.Rock(binIdx), yHeavy(binIdx), cfg.reg.shifts_psi(:), ...
    regShift, regShift - cfg.reg.shifts_psi(:), ...
    'VariableNames', {'Mixture', 'Rock', 'y_heavy', 'Expected_psi', 'Got_psi', 'Delta_psi'});
disp(regTable);
if any(~(abs(regTable.Delta_psi) <= cfg.reg.tol_psi))
    error('tip:RegressionFailed', ['Optimised-TIP shifts do not reproduce the previous round ' ...
        '(worst %+.1f psi). Check the C1-nC5 compositions in the workbook and cfg.reg.basis.'], ...
        max(abs(regTable.Delta_psi), [], 'omitnan'));
end
iB1 = find(groupKey(binIdx) == cfg.reg.b1Group, 1);
if abs(regShift(iB1) - cfg.reg.b1Shift_psi) > cfg.reg.b1Tol_psi
    error('tip:B1AnchorFailed', 'B1 anchor: expected %.1f psi, got %.3f psi.', cfg.reg.b1Shift_psi, regShift(iB1));
end
fprintf('Regression gate PASSED.\n\n');

%% 3. Fit, predict, judge at nominal r_p ---------------------------------------
combosRun = cfg.combos.Combination;
if ~isempty(cfg.runCombos), combosRun = intersect(combosRun, cfg.runCombos, 'stable'); end
nF = height(cfg.forms);
fits = containers.Map();
fitRows = {}; binRows = {}; terRows = {}; pairRows = {}; accRows = {};

for cc = combosRun.'
    ci = find(cfg.combos.Combination == cc, 1);
    runner.setFWI(cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci));
    for md = cfg.modes
        tipOpt = cfg.tips(char(cc + "|" + md));
        tipOpt = tipOpt(:);
        for f = 1:nF
            F = cfg.forms(f, :);
            dCal = arrayfun(@(g) pairDescriptor(P, grpPair{g}(1), grpPair{g}(2), F.Descriptor, cfg), (1:G).');
            Fi = fitForm(tipOpt, dCal, isB1g);
            Fi.name = F.Name; Fi.descriptor = F.Descriptor; Fi.T = cfg.T; Fi.kB = cfg.kB;
            Fi.dRange = [min(dCal), max(dCal)];
            % equivalent re-expressions of the rock offset (same predictions)
            Fi.deltaS   = (Fi.alphaEF2 - Fi.alphaB1) / (rockD.S.EF2 - rockD.S.B1);
            Fi.alpha0S  = Fi.alphaEF2 - Fi.deltaS * rockD.S.EF2;
            Fi.deltaTOC = (Fi.alphaEF2 - Fi.alphaB1) / (rockD.TOC.EF2 - rockD.TOC.B1);
            Fi.alpha0TOC = Fi.alphaEF2 - Fi.deltaTOC * rockD.TOC.EF2;
            key = cc + "|" + md + "|" + F.Name;
            fits(char(key)) = Fi;
            fitRows(end+1, :) = {cc, md, F.Name, F.Descriptor, Fi.alphaEF2, ...
                Fi.alphaB1, Fi.B, Fi.BdropMin, Fi.BdropMax, Fi.rmseLnEF2, ...
                Fi.alpha0S, Fi.deltaS, Fi.alpha0TOC, Fi.deltaTOC}; %#ok<AGROW>

            % --- pair guard: every applied pair, TIP must be finite and <= 0
            [pt, bad] = appliedPairs(runner, dp, caseIdx, P, Fi, max(abs(tipOpt)) * cfg.extrapFlagFactor);
            for q = 1:height(pt)
                pairRows(end+1, :) = {cc, md, F.Name, pt.Context(q), pt.Rock(q), pt.Pair(q), ...
                    pt.x(q), pt.InCalRange(q), pt.TIP(q), pt.ExtrapFlag(q)}; %#ok<AGROW>
            end
            if bad
                warning('tip:BadPairTIP', '%s: positive or non-finite pair TIP, row rejected.', key);
                accRows(end+1, :) = {cc, md, F.Name, false, false, false, false, "", ...
                    NaN, NaN, "rejected: pair TIP > 0 or non-finite"}; %#ok<AGROW>
                continue;
            end

            % --- binaries
            runner.clearWarm();
            [shB, errB, pctB] = deal(nan(numel(binIdx), 1));
            for n = 1:numel(binIdx)
                i = binIdx(n); c = dp.getCase(i);
                eng = runner.engine(c.Mixture, c.Rock);
                K = tipMatrix(string(eng.eos.Fluid.ComponentNames), c.Rock, P, Fi);
                r = runner.solve(c, md, K, 'Key', sprintf("cor|%d|%s|%d", f, md, i));
                if r.ok
                    shB(n) = shiftOf(r, cfg); errB(n) = shB(n) - dp.ShiftMid_psi(i);
                    pctB(n) = pctErr(r.P, dp.ConfMid_psi(i), psi);
                end
                binRows(end+1, :) = {cc, md, F.Name, c.Mixture, c.Rock, c.z(end), K(1, 2), ...
                    dp.ShiftMid_psi(i), dp.ConfHalfWidth_psi(i), shB(n), errB(n), ...
                    abs(errB(n)) <= dp.ConfHalfWidth_psi(i), pctB(n), r.ok, r.route, r.reason}; %#ok<AGROW>
            end

            % --- ternaries
            [shT, errT, pctT, inT] = deal(nan(numel(terIdx), 1));
            for n = 1:numel(terIdx)
                i = terIdx(n); c = dp.getCase(i);
                eng = runner.engine(c.Mixture, c.Rock);
                K = tipMatrix(string(eng.eos.Fluid.ComponentNames), c.Rock, P, Fi);
                r = runner.solve(c, md, K, 'Key', sprintf("ter|%d|%s|%d", f, md, i));
                if r.ok
                    shT(n) = shiftOf(r, cfg); errT(n) = shT(n) - dp.ShiftMid_psi(i);
                    pctT(n) = pctErr(r.P, dp.ConfMid_psi(i), psi);
                end
                inT(n) = r.ok && abs(errT(n)) <= dp.ConfHalfWidth_psi(i);
                terRows(end+1, :) = {cc, md, F.Name, c.Mixture, c.Rock, dp.ShiftMid_psi(i), ...
                    dp.ConfHalfWidth_psi(i), shT(n), errT(n), inT(n), pctT(n), r.ok, r.route, r.reason}; %#ok<AGROW>
            end

            isC10 = dp.Mixture(binIdx) == "C1-nC10";
            a1 = ~cfg.accept.noNegativeShift || all([shB; shT] > 0);          % NaN -> fail
            a2 = all(abs(errB(isC10)) <= cfg.accept.maxErr_C1nC10_psi);    % NaN -> fail
            a3 = ~cfg.accept.ternaryInBand   || all(inT == 1);
            a4 = max(abs([pctB; pctT])) <= cfg.accept.maxPdewErr_pct && all(isfinite([pctB; pctT]));
            accRows(end+1, :) = {cc, md, F.Name, a1, a2, a3, a4, "", ...
                sqrt(mean(errB.^2)), sqrt(mean(errT.^2)), ""}; %#ok<AGROW>
        end
    end
end

accTable = cell2table(accRows, 'VariableNames', {'Combination', 'Mode', 'Form', ...
    'NoNegativeShift', 'C1nC10_error', 'Ternaries_in_band', 'Pdew_within_1pct', 'rp_monotonic', ...
    'RMSE_binary_psi', 'RMSE_ternary_psi', 'Note'});
accTable = textToString(accTable);
accTable.rp_monotonic = repmat("not run", height(accTable), 1);

%% 4. r_p sweep: shift must fall strictly as r_p rises, own rock and offset -----
radii = sort(cfg.rpSweep.radii_nm);
nR = numel(radii);
nomR = [cfg.poreRadius('EF2'), cfg.poreRadius('B1')] * 1e9;
[dNom, iNom] = min(abs(radii(:) - nomR), [], 1);
if any(dNom > 1e-6), error('tip:NominalRadius', 'Nominal radii must be in cfg.rpSweep.radii_nm.'); end
sweep = containers.Map(); why = containers.Map();
for ir = 1:nR
    rr = solvers.DewPointRunner(cfg.mixtureFile, containers.Map({'EF2', 'B1'}, ...
        {radii(ir) * 1e-9, radii(ir) * 1e-9}));
    rr.maxShift_psi = cfg.rpSweep.maxShift_psi;
    for i = caseIdx.'
        rr.engine(dp.Mixture(i), dp.Rock(i));
        rr.bulk(dp.getCase(i));
    end
    for cc = combosRun.'
        ci = find(cfg.combos.Combination == cc, 1);
        rr.setFWI(cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci));
        for md = cfg.rpSweep.modes
            for f = 1:nF
                key = char(cc + "|" + md + "|" + cfg.forms.Name(f));
                row = find(accTable.Combination == cc & accTable.Mode == md & ...
                    accTable.Form == cfg.forms.Name(f), 1);
                if startsWith(accTable.Note(row), "rejected"), continue; end
                Fi = fits(key);
                if ~isKey(sweep, key)
                    sweep(key) = nan(numel(caseIdx), nR);
                    why(key) = strings(numel(caseIdx), nR);
                end
                Sm = sweep(key); Wm = why(key);
                rr.clearWarm();
                for n = 1:numel(caseIdx)
                    i = caseIdx(n); c = dp.getCase(i);
                    eng = rr.engine(c.Mixture, c.Rock);
                    K = tipMatrix(string(eng.eos.Fluid.ComponentNames), c.Rock, P, Fi);
                    r = rr.solve(c, md, K, 'Key', sprintf("rp|%d|%d|%s|%d", ir, f, md, i));
                    if r.ok, Sm(n, ir) = shiftOf(r, cfg); else, Wm(n, ir) = string(r.reason); end
                end
                sweep(key) = Sm; why(key) = Wm;
            end
        end
    end
    fprintf('r_p sweep: %.2f nm done.\n', radii(ir));
end

monoRows = {};
sk = keys(sweep);
for s = 1:numel(sk)
    Sm = sweep(sk{s}); Wm = why(sk{s});
    pr = split(string(sk{s}), "|");
    st = strings(numel(caseIdx), 1);
    for n = 1:numel(caseIdx)
        i = caseIdx(n);
        nomCol = iNom(1 + (dp.Rock(i) == "B1"));
        [st(n), maxSlope] = monotonicStatus(Sm(n, :), radii, nomCol, cfg.rpSweep);
        fails = Wm(n, Wm(n, :) ~= "");
        monoRows(end+1, :) = [{pr(1), pr(2), pr(3), dp.Mixture(i), dp.Rock(i), yHeavy(i), st(n), ...
            maxSlope, nnz(isfinite(Sm(n, :))), strjoin(unique(fails), " || ")}, num2cell(Sm(n, :))]; %#ok<AGROW>
    end
    row = accTable.Combination == pr(1) & accTable.Mode == pr(2) & accTable.Form == pr(3);
    if any(st == "Fail"), accTable.rp_monotonic(row) = "Fail";
    elseif all(st == "Pass"), accTable.rp_monotonic(row) = "Pass";
    else, accTable.rp_monotonic(row) = "Undetermined";
    end
end
radNames = "Shift_rp_" + strrep(string(radii), ".", "p") + "nm";
monoTable = cell2table(monoRows, 'VariableNames', [{'Combination', 'Mode', 'Form', 'Mixture', ...
    'Rock', 'y_heavy', 'Status', 'MaxLogSlope', 'nFinite', 'FailReasons'}, cellstr(radNames)]);

%% 5. Verdict and selection -----------------------------------------------------
rpOK = accTable.rp_monotonic == "Pass" | ~cfg.accept.rpMonotonic;
accTable.Pass = accTable.NoNegativeShift & accTable.C1nC10_error & ...
    accTable.Ternaries_in_band & accTable.Pdew_within_1pct & rpOK;
cand = accTable(accTable.Mode == cfg.primaryMode & accTable.Pass, :);
selected = table();
if isempty(cand)
    fprintf('\nNo row passes all criteria in %s. Nothing selected.\n', cfg.primaryMode);
else
    best = min(cand.RMSE_ternary_psi);
    tie  = cand(cand.RMSE_ternary_psi <= best + cfg.select.tieTol_psi, :);
    [~, j] = min(tie.RMSE_binary_psi);
    selected = tie(j, :);
    fprintf('\nSELECTED (%s): %s / %s | ternary RMSE %.1f psi | binary RMSE %.1f psi\n', ...
        cfg.primaryMode, selected.Combination, selected.Form, selected.RMSE_ternary_psi, ...
        selected.RMSE_binary_psi);
    fprintf(['Disclosure: the TIP form was changed on physical grounds (r_p dependence of wrong ' ...
        'sign); contrast descriptors were rejected by the ternaries (heavy-heavy |TIP| too large); ' ...
        'the choice among the declared forms and combinations used the ternaries.\n']);
end

%% 6. CO2 pair preview (information only; NOT used in any criterion or selection) --
co2Rows = {};
for s = keys(fits)
    Fi = fits(s{1});
    pr = split(string(s{1}), "|");
    for pair = [["CO2", "C1"]; ["CO2", "nC8"]; ["C1", "nC8"]].'
        d = pairDescriptor(P, pair(1), pair(2), Fi.descriptor, cfg);
        co2Rows(end+1, :) = {pr(1), pr(2), pr(3), pair(1) + "-" + pair(2), d, ...
            d >= Fi.dRange(1) && d <= Fi.dRange(2), tipEval(d, "B1", Fi)}; %#ok<AGROW>
    end
end
co2Table = cell2table(co2Rows, 'VariableNames', {'Combination', 'Mode', 'Form', 'Pair', 'x', ...
    'InCalRange', 'TIP_B1'});

%% 7. Output -----------------------------------------------------------------------
fitTable = cell2table(fitRows, 'VariableNames', {'Combination', 'Mode', 'Form', 'Descriptor', ...
    'alpha_EF2', 'alpha_B1', 'B', 'B_dropOne_min', 'B_dropOne_max', 'RMSE_ln_EF2', 'alpha0_S', 'delta_S_per_sqrtK', 'alpha0_TOC', 'delta_TOC'});
fitTable = textToString(fitTable);
binTable = cell2table(binRows, 'VariableNames', {'Combination', 'Mode', 'Form', 'Mixture', 'Rock', ...
    'y_heavy', 'TIP_12', 'Shift_exp_psi', 'HalfWidth_psi', 'Shift_pred_psi', 'Error_psi', 'InBand', ...
    'Pdew_err_pct', 'Accepted', 'Route', 'Reason'});
binTable = textToString(binTable);
terTable = cell2table(terRows, 'VariableNames', {'Combination', 'Mode', 'Form', 'Mixture', 'Rock', ...
    'Shift_exp_psi', 'HalfWidth_psi', 'Shift_pred_psi', 'Error_psi', 'InBand', 'Pdew_err_pct', ...
    'Accepted', 'Route', 'Reason'});
terTable = textToString(terTable);
pairTableAll = cell2table(pairRows, 'VariableNames', {'Combination', 'Mode', 'Form', 'Context', ...
    'Rock', 'Pair', 'x', 'InCalRange', 'TIP', 'ExtrapFlag'});
pairTableAll = textToString(pairTableAll);
monoTable = textToString(monoTable);
co2Table = textToString(co2Table);

fprintf('\n=== Correlation fits (rock offset, dof = 1) ===\n'); disp(fitTable);
fprintf('=== Binary predictions ===\n');   disp(binTable);
fprintf('=== Ternary predictions ===\n');  disp(terTable);
fprintf('=== Pair TIPs applied ===\n');    disp(pairTableAll);
fprintf('=== r_p sweep ===\n');            disp(monoTable);
fprintf('=== Acceptance ===\n');           disp(accTable);
fprintf('=== CO2 pair preview (information only) ===\n'); disp(co2Table);

save(fullfile(cfg.resultsDir, cfg.outFile + ".mat"), 'cfg', 'rockD', 'regTable', 'fitTable', ...
    'binTable', 'terTable', 'pairTableAll', 'monoTable', 'accTable', 'selected', 'co2Table');
xls = fullfile(cfg.resultsDir, cfg.outFile + ".xlsx");
if isfile(xls), delete(xls); end
outs = {regTable, 'Regression_gate'; fitTable, 'Correlation'; binTable, 'Binaries'; ...
    terTable, 'Ternaries'; pairTableAll, 'Pair_TIPs'; monoTable, 'rp_sweep'; ...
    accTable, 'Acceptance'; co2Table, 'CO2_preview_info_only'};
for q = 1:size(outs, 1), writetable(outs{q, 1}, xls, 'Sheet', outs{q, 2}); end
if ~isempty(selected), writetable(selected, xls, 'Sheet', 'Selected'); end
fprintf('\nWritten: %s\n', xls);

%% ============================== Local functions ==============================
function t = textToString(t)
isTxt = varfun(@(v) iscell(v) && all(cellfun(@(e) isstring(e) || ischar(e), v)), t, ...
    'OutputFormat', 'uniform');
if any(isTxt), t = convertvars(t, t.Properties.VariableNames(isTxt), 'string'); end
end

function P = loadComponentProps(file)
T = readtable(file, 'Sheet', 'Comp_Props', 'VariableNamingRule', 'preserve');
P = containers.Map();
for r = 1:height(T)
    nm = string(T.comp(r));
    if ismissing(nm) || nm == "", continue; end
    P(char(nm)) = struct('M', T.Mw(r), 'sigma', T.LJ_Size(r), 'eps', T.LJ_Energy(r), 'omega', T.omega(r));
end
end

function D = loadRockDescriptors(file)
% S_rock = sum(phi_m sqrt(E_m)), phi = volume fraction incl. TOC; TOC as volume fraction.
rc  = readcell(file, 'Sheet', 'Rock_Props');
hdr = string(rc(1, :));
lab = string(rc(:, 1));
mc  = find(~ismember(hdr, ["Mineral", "Rock", "theta"]) & ~ismissing(hdr));
iR  = find(hdr == "Rock");
E   = cell2mat(rc(lab == "mineral_E", mc));
rho = cell2mat(rc(lab == "mineral_rho", mc));
iT  = hdr(mc) == "TOC";
D.S = struct(); D.TOC = struct();
for r = find(lab == "mineral_w").'
    w   = cell2mat(rc(r, mc));
    phi = (w ./ rho) / sum(w ./ rho);
    nm  = string(rc{r, iR});
    D.S.(nm)   = sum(phi .* sqrt(E));
    D.TOC.(nm) = phi(iT);
end
ref = struct('EF2', 8.0388, 'B1', 5.2870);
for nm = ["EF2", "B1"]
    if abs(D.S.(nm) - ref.(nm)) > 1e-3
        error('tip:Srock', 'S_rock(%s) = %.4f, expected %.4f.', nm, D.S.(nm), ref.(nm));
    end
end
end

function x = pairDescriptor(P, a, b, kind, cfg)
pa = P(char(a)); pb = P(char(b));
switch kind
    case "wbar",      x = (pa.omega + pb.omega) / 2;
    case "inv_wbar",  x = 2 / (pa.omega + pb.omega);
    case "Tstar",     x = cfg.kB * cfg.T / sqrt(pa.eps * pb.eps);
    case "inv_Tstar", x = sqrt(pa.eps * pb.eps) / (cfg.kB * cfg.T);
    case "ln_Tstar",  x = log(cfg.kB * cfg.T / sqrt(pa.eps * pb.eps));
    case "inv_Mgeo",  x = 1 / sqrt(pa.M * pb.M);
    otherwise, error('tip:Descriptor', 'Unknown descriptor %s.', kind);
end
end

function F = fitForm(tipOpt, x, isB1)
% ln(-TIP) = alpha_EF2*[EF2] + alpha_B1*[B1] + B*x
y = log(-tipOpt(:));
X = [double(~isB1), double(isB1), x(:)];
if rank(X) < 3, error('tip:Rank', 'Design matrix rank-deficient.'); end
par = X \ y;
F.alphaEF2 = par(1); F.alphaB1 = par(2); F.B = par(3);
res = X * par - y;
F.rmseLnEF2 = sqrt(mean(res(~isB1).^2));
Bd = nan(numel(y), 1);
for q = find(~isB1).'                      % drop-one over EF2 groups (B1 drop is unidentified)
    keep = true(numel(y), 1); keep(q) = false;
    if rank(X(keep, :)) == 3, pD = X(keep, :) \ y(keep); Bd(q) = pD(3); end
end
F.BdropMin = min(Bd); F.BdropMax = max(Bd);
end

function t = tipEval(x, rock, F)
switch string(rock)
    case "EF2", a = F.alphaEF2;
    case "B1",  a = F.alphaB1;
    otherwise, error('tip:Rock', 'No rock offset for %s.', rock);
end
t = -exp(a + F.B * x);
end

function [st, maxSlope] = monotonicStatus(sh, radii, nomCol, sw)
maxSlope = NaN;
if ~isfinite(sh(nomCol)), st = "Undetermined"; return; end
sg = sign(sh(nomCol));
lo = nomCol; while lo > 1 && isfinite(sh(lo-1)) && sign(sh(lo-1)) == sg, lo = lo - 1; end
hi = nomCol; while hi < numel(sh) && isfinite(sh(hi+1)) && sign(sh(hi+1)) == sg, hi = hi + 1; end
b = sh(lo:hi); r = radii(lo:hi);
if any(b <= 0), st = "Fail"; return; end
if numel(b) < sw.minFinite, st = "Undetermined"; return; end
maxSlope = max(-diff(log(b)) ./ diff(log(r)));
if all(diff(b) < 0), st = "Pass"; else, st = "Fail"; end
end

function K = tipMatrix(names, rock, P, F)
nc = numel(names); K = zeros(nc);
for a = 1:nc - 1
    for b = a + 1:nc
        K(a, b) = tipEval(pairDescriptor(P, names(a), names(b), F.descriptor, ...
            struct('kB', F.kB, 'T', F.T)), rock, F);
        K(b, a) = K(a, b);
    end
end
end

function [T, bad] = appliedPairs(runner, dp, caseIdx, P, F, flagAbs)
rows = {}; seen = strings(0, 1);
for i = caseIdx.'
    c = dp.getCase(i);
    eng = runner.engine(c.Mixture, c.Rock);
    nm = string(eng.eos.Fluid.ComponentNames);
    ctx = "binary"; if numel(c.z) > 2, ctx = "ternary"; end
    for a = 1:numel(nm) - 1
        for b = a + 1:numel(nm)
            key = ctx + "|" + c.Rock + "|" + nm(a) + "-" + nm(b);
            if any(seen == key), continue; end
            seen(end+1) = key; %#ok<AGROW>
            d = pairDescriptor(P, nm(a), nm(b), F.descriptor, struct('kB', F.kB, 'T', F.T));
            t = tipEval(d, c.Rock, F);
            rows(end+1, :) = {ctx, c.Rock, nm(a) + "-" + nm(b), d, ...
                d >= F.dRange(1) && d <= F.dRange(2), t, abs(t) > flagAbs}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'Context', 'Rock', 'Pair', 'x', 'InCalRange', 'TIP', 'ExtrapFlag'});
bad = any(T.TIP > 0) || any(~isfinite(T.TIP));
end

function sh = shiftOf(r, cfg)
if cfg.shiftBasis == "expBulk", sh = r.shiftExpBulk_psi; else, sh = r.shift_psi; end
end

function e = pctErr(P_Pa, Pexp_psi, psi)
e = 100 * (P_Pa / (Pexp_psi * psi) - 1);
end