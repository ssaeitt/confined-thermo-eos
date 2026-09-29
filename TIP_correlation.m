%% TIP_CORRELATION  ln(-TIP_ij) = alpha + beta / M_ij + gamma * (sigma_bar_ij / r_p)
%
%   TIP_ij = -exp( alpha + beta / M_ij + gamma * sigma_bar_ij / r_p )
%
% Sign-safe by construction: TIP < 0 for every pair, including pairs the
% correlation is only applied to (ternary pairs, later CO2 pairs).
% M_ij: "heavy" = max(M_i, M_j) is primary; "geomean" = sqrt(M_i M_j) is the
% sensitivity case. A linear form is kept behind cfg.tipForm for comparison.
%
% Three parameters on four binary groups: ONE degree of freedom. gamma is a
% TWO-ROCK INTERPOLATION - only EF2 (54 nm) and B1 (11.25 nm) are sampled, so
% it is a contrast between two pore sizes, not a trend in 1/r_p.
%
% WALLS: each mixture engine keeps its own rock mineralogy (EF2 / B1). The
% quartz-only wall belongs to the dTc calibration in Fit_parameters.m and is
% NOT applied here.
%
% Nothing is reported until the regression gate in section 2 passes.

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
cfg.mixtureFile = fullfile(scriptDir, "config", "MixtureData.xlsx");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.outFile     = "tip_correlation";
cfg.T           = 293.15;
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});   % [m], radii
cfg.modes       = ["FWIOnly", "Combined"];
cfg.tipForm     = "loglinear";          % "loglinear" | "linear"
cfg.shiftBasis  = "modelBulk";          % basis for reported/compared shifts:
                                        % "modelBulk" removes the bulk-EOS bias,
                                        % "expBulk" reproduces the V31 definition
cfg.Mdefs       = ["heavy", "geomean"]; % first entry is primary

% Locked (k, p, lambda) from Fit_parameters; sequential-exact first.
% NOTE: the reported joint-exact k was 4.51381302592871, which is the
% sequential-exact LAMBDA - a transcription slip. Its p and lambda match the
% earlier joint-exact run exactly, so k is taken from that run (98.05788478).
% Correct it here if the true value differs.
cfg.combos = table( ...
    ["sequential-exact"; "joint-exact"; "joint-linearized"], ...
    [151.385346647943;  98.0578847776961;  112.42634826677 ], ...
    [1.30682834548441;  1.69745310322519;  1.3558179290678 ], ...
    [4.51381302592871;  10.2425163919516;  2.0916235090532 ], ...
    'VariableNames', {'Combination', 'k', 'p', 'lambda'});
cfg.runCombos = [];                     % subset to run; [] or "" = all three

% Sanity check on the locked parameters (catches transcription slips)
bad = cfg.combos.k < 10 | cfg.combos.k > 1e4 | cfg.combos.p < 0.5 | cfg.combos.p > 3 | ...
      cfg.combos.lambda < 0 | cfg.combos.lambda > 50;
if any(bad)
    error('tip:ImplausibleParameters', 'Check (k, p, lambda) for: %s', ...
          strjoin(cfg.combos.Combination(bad), ", "));
end

cfg.groups = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"];
cfg.tips = containers.Map( ...
    {'sequential-exact|FWIOnly',  'joint-exact|FWIOnly',  'joint-linearized|FWIOnly', ...
     'sequential-exact|Combined', 'joint-exact|Combined', 'joint-linearized|Combined'}, ...
    {[-3.48386287716462,  -1.46205964046861, -1.12088738883485, -0.905508761783693], ...
     [-13.4437287153297,  -5.90838752675671, -4.59334889956389, -2.38867380361359 ], ...
     [-4.9229696302796,   -1.96920366521686, -1.44628069292089, -0.754862049795788], ...
     [-3.40199024752739,  -1.45777900094979, -1.11884643084433, -0.899502001873144], ...
     [-14.4508613237875,  -5.84525317456513, -4.79513444682033, -2.24015630615276 ], ...
     [-4.76033469976459,  -1.96111617688159, -1.44196236936452, -0.74629351212583 ]});

% --- Regression gate: must pass before any result is read -------------------
cfg.reg.combination = "sequential-exact";
cfg.reg.mode        = "FWIOnly";
cfg.reg.shifts_psi  = [110.9, 39.1, -8.9, 175.0, 140.5, 84.4, 98.0, 313.7, 180.6];  % Mix_Props order
cfg.reg.tol_psi     = 1.0;
cfg.reg.basis       = "expBulk";        % basis of cfg.reg.shifts_psi:
                                        %   "expBulk"   P_model - P_bulk,EXPERIMENTAL (V31 basis,
                                        %               used by Fit_parameters' Shift_model_psi)
                                        %   "modelBulk" P_model - P_bulk,MODEL
cfg.reg.b1Group     = "C1-nC8|B1";
cfg.reg.b1Shift_psi = 98.0;
cfg.reg.b1Tol_psi   = 0.1;              % "exactly": optimised TIP reproduces the B1 anchor

% --- Acceptance criteria, declared BEFORE the run ---------------------------
cfg.accept.noNegativeC1nC5   = true;    % no negative predicted C1-nC5 shift
cfg.accept.maxErr_C1nC10_psi = 60;      % C1-nC10 over-prediction must drop from ~+170 psi
cfg.accept.ternaryInBand     = true;    % both ternaries inside their experimental band
cfg.accept.maxPdewErr_pct    = 1.0;     % absolute dew-point error (supervisor's target)

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
psi = entities.DewPointDataset.psi2Pa;

fprintf('Acceptance criteria (declared before the run):\n');
fprintf('  no negative C1-nC5 shift: %d | max C1-nC10 error <= %.0f psi | both ternaries in band: %d | |Pdew err| <= %.1f %%\n', ...
        cfg.accept.noNegativeC1nC5, cfg.accept.maxErr_C1nC10_psi, cfg.accept.ternaryInBand, ...
        cfg.accept.maxPdewErr_pct);

%% 1. Data, engines (own mineralogy), regressors -----------------------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius);

groupKey = dp.Mixture + "|" + dp.Rock;
isBin  = dp.HasConfined & dp.NumComponents == 2 & ismember(groupKey, cfg.groups);
isTer  = dp.HasConfined & dp.NumComponents > 2;
binIdx = find(isBin);
for i = [binIdx; find(isTer)].'
    runner.engine(dp.Mixture(i), dp.Rock(i));
end

% Wall check: every engine must carry its own rock, not a calibration wall
ek = keys(runner.Engines);
for e = 1:numel(ek)
    eng = runner.Engines(ek{e});
    pr  = split(string(ek{e}), "|");
    rk  = eng.eos.Rock;
    if rk.RockName ~= pr(2) || rk.NumMinerals < 2
        error('tip:WallMismatch', ['Engine %s carries wall "%s" with %d mineral(s). The ' ...
            'quartz-only wall is for the dTc calibration only; mixture engines must use ' ...
            'the full mineralogy of their rock.'], string(ek{e}), rk.RockName, rk.NumMinerals);
    end
    fprintf('Wall check: %-14s -> %s, %d minerals, TOC %.4f\n', string(ek{e}), rk.RockName, ...
            rk.NumMinerals, rk.RawTOC);
end

for i = [binIdx; find(isTer)].'
    [~, ~, okb] = runner.bulk(dp.getCase(i));
    if ~okb
        error('tip:BulkFailed', 'Bulk dew point failed for %s (%s).', dp.Mixture(i), dp.Rock(i));
    end
end

G = numel(cfg.groups);
[sigOverR, Mheavy, Mgeo] = deal(zeros(G, 1));
for g = 1:G
    pr  = split(cfg.groups(g), "|");
    eng = runner.engine(pr(1), pr(2));
    sig = eng.eos.Fluid.LJ_Size(:);
    MW  = eng.eos.Fluid.MW(:);
    sigOverR(g) = 0.5 * (sig(1) + sig(2)) / cfg.poreRadius(char(pr(2)));
    Mheavy(g)   = max(MW);
    Mgeo(g)     = sqrt(prod(MW));
end
fprintf('\nGroup regressors:\n');
disp(table(cfg.groups.', sigOverR, Mheavy, Mgeo, ...
     'VariableNames', {'Group', 'sigma_over_rp', 'M_heavy', 'M_geomean'}));

%% 2. Regression gate ---------------------------------------------------------
fprintf('=== Regression gate (%s / %s, optimised TIPs) ===\n', cfg.reg.combination, cfg.reg.mode);
ci = find(cfg.combos.Combination == cfg.reg.combination, 1);
runner.setFWI(cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci));
runner.clearWarm();
tipReg = cfg.tips(char(cfg.reg.combination + "|" + cfg.reg.mode)).';
[regModel, regExp] = deal(nan(numel(binIdx), 1));
for n = 1:numel(binIdx)
    i = binIdx(n);
    c = dp.getCase(i);
    g = find(cfg.groups == groupKey(i), 1);
    r = runner.solve(c, cfg.reg.mode, [0, tipReg(g); tipReg(g), 0], 'Key', sprintf("reg|%d", i));
    if r.ok
        regModel(n) = r.shift_psi;           % vs model bulk
        regExp(n)   = r.shiftExpBulk_psi;    % vs experimental bulk (V31)
    end
end
regShift = regExp;
if cfg.reg.basis == "modelBulk", regShift = regModel; end
regTable = table(dp.Mixture(binIdx), dp.Rock(binIdx), cellfun(@(v) v(end), dp.MoleFrac(binIdx)), ...
    cfg.reg.shifts_psi(:), regShift, regShift - cfg.reg.shifts_psi(:), regExp, regModel, ...
    regModel - regExp, ...
    'VariableNames', {'Mixture', 'Rock', 'y_heavy', 'Expected_psi', 'Got_psi', 'Delta_psi', ...
    'Shift_vs_expBulk', 'Shift_vs_modelBulk', 'BulkBias_psi'});
disp(regTable);
bad = ~(abs(regTable.Delta_psi) <= cfg.reg.tol_psi);
if any(bad)
    other = regModel;
    if cfg.reg.basis == "modelBulk", other = regExp; end
    dOther = max(abs(other - cfg.reg.shifts_psi(:)), [], 'omitnan');
    hint = "";
    if dOther <= cfg.reg.tol_psi
        hint = sprintf([' The OTHER basis matches within %.1f psi: set cfg.reg.basis to the one ' ...
            'the expected numbers were produced with.'], cfg.reg.tol_psi);
    else
        hint = sprintf(' The other basis is off by up to %+.1f psi as well.', dOther);
    end
    error('tip:RegressionFailed', ['Optimised-TIP shifts do not reproduce the previous round ' ...
        'on basis "%s" (%d of %d rows off by more than %.1f psi, worst %+.1f).%s Results are ' ...
        'not readable until this passes.'], cfg.reg.basis, nnz(bad), height(regTable), ...
        cfg.reg.tol_psi, max(abs(regTable.Delta_psi), [], 'omitnan'), hint);
end
iB1 = find(groupKey(binIdx) == cfg.reg.b1Group);
if isempty(iB1) || abs(regShift(iB1(1)) - cfg.reg.b1Shift_psi) > cfg.reg.b1Tol_psi
    error('tip:B1AnchorFailed', 'B1 anchor: expected %.1f psi, got %.3f psi.', ...
          cfg.reg.b1Shift_psi, regShift(iB1(1)));
end
fprintf('Regression gate PASSED (B1 anchor %.3f psi, worst delta %+.2f psi).\n\n', ...
        regShift(iB1(1)), max(abs(regTable.Delta_psi), [], 'omitnan'));

%% 3. Fit, predict, judge ------------------------------------------------------
combosRun = cfg.combos.Combination;
if ~isempty(cfg.runCombos) && cfg.runCombos ~= ""
    combosRun = intersect(combosRun, cfg.runCombos, 'stable');
end

fitRows = {};  levRows = {};  binRows = {};  terRows = {};  gapRows = {};  pairRows = {};  accRows = {};
for cc = combosRun.'
    ci = find(cfg.combos.Combination == cc, 1);
    runner.setFWI(cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci));

    for md = cfg.modes
        tipOpt = cfg.tips(char(cc + "|" + md)).';

        % --- calibration baseline: per-group optimised TIPs
        runner.clearWarm();
        errOpt = nan(numel(binIdx), 1);  okOpt = false(numel(binIdx), 1);
        for n = 1:numel(binIdx)
            i = binIdx(n);
            c = dp.getCase(i);
            g = find(cfg.groups == groupKey(i), 1);
            r = runner.solve(c, md, [0, tipOpt(g); tipOpt(g), 0], 'Key', sprintf("opt|%d|%s", i, md));
            sh = shiftOf(r, cfg);
            errOpt(n) = sh - dp.ShiftMid_psi(i);
            okOpt(n)  = r.ok;
            binRows(end+1, :) = {cc, md, "optimised", c.Mixture, c.Rock, c.z(end), tipOpt(g), ...
                tipOpt(g), dp.ShiftMid_psi(i), dp.ConfHalfWidth_psi(i), sh, errOpt(n), ...
                abs(errOpt(n)) <= dp.ConfHalfWidth_psi(i), pctErr(r.P, dp.ConfMid_psi(i), psi), ...
                r.ok, r.route, r.reason}; %#ok<AGROW>
        end
        rmseOpt = sqrt(mean(errOpt(okOpt).^2));

        for dI = 1:numel(cfg.Mdefs)
            Mdef = cfg.Mdefs(dI);
            role = "sensitivity";
            if dI == 1, role = "primary"; end
            Mg = Mheavy;
            if Mdef == "geomean", Mg = Mgeo; end

            % --- fit: ln(-TIP) = alpha + beta/M + gamma (sigma_bar/rp)
            X = [ones(G, 1), 1 ./ Mg, sigOverR];
            if cfg.tipForm == "loglinear"
                if any(tipOpt >= 0)
                    error('tip:PositiveTIP', 'Log-linear form needs TIP < 0 for all fitted groups.');
                end
                yFit = log(-tipOpt);          % fit is in ln(-TIP); TIP = -exp(...) is sign-safe
            else
                yFit = tipOpt;
            end
            par   = X \ yFit;
            tipOf = makeTipFun(par, cfg.tipForm);
            resid = X * par - yFit;
            dof = G - 3;                                          % = 1
            s2  = sum(resid.^2) / max(dof, 1);
            C   = s2 * inv(X.' * X);                              %#ok<MINV>
            se  = sqrt(diag(C));
            lev = diag(X * ((X.' * X) \ X.'));
            isB1g = endsWith(cfg.groups(:), "|B1");
            tipCor = arrayfun(@(g) tipOf(Mg(g), sigOverR(g)), (1:G).');
            rmseTip = sqrt(mean((tipCor - tipOpt).^2));

            for g = 1:G
                keep = true(G, 1);  keep(g) = false;
                gDrop = NaN;
                if rank(X(keep, :)) == 3
                    pD = X(keep, :) \ yFit(keep);
                    gDrop = pD(3);
                end
                levRows(end+1, :) = {cc, md, Mdef, role, cfg.groups(g), sigOverR(g), Mg(g), lev(g), ...
                    tipOpt(g), tipCor(g), tipCor(g) - tipOpt(g), gDrop, gDrop - par(3)}; %#ok<AGROW>
            end
            fitRows(end+1, :) = {cc, md, Mdef, role, cfg.tipForm, par(1), se(1), par(2), se(2), ...
                par(3), se(3), rmseTip, sqrt(mean(resid.^2)), dof, max(lev(isB1g)), ...
                "gamma = two-rock interpolation (EF2 54 nm, B1 11.25 nm), dof = 1"}; %#ok<AGROW>
            fprintf(['\n%s / %s / M = %s (%s):  alpha = %+.4f (se %.4f), beta = %+.4f (se %.4f), ' ...
                     'gamma = %+.4f (se %.4f)\n  RMSE_TIP = %.4f, dof = %d, B1 leverage = %.3f ' ...
                     '| gamma: two-rock interpolation (EF2, B1)\n'], cc, md, Mdef, role, par(1), ...
                    se(1), par(2), se(2), par(3), se(3), rmseTip, dof, max(lev(isB1g)));

            % --- pair-TIP guard: print every pair used, reject TIP >= 0
            [pr, bad] = pairTable(runner, dp, cfg, binIdx, isTer, Mdef, tipOf);
            for q = 1:height(pr)
                pairRows(end+1, :) = {cc, md, Mdef, pr.Context(q), pr.Rock(q), pr.Pair(q), ...
                    pr.M(q), pr.sigma_over_rp(q), pr.TIP(q)}; %#ok<AGROW>
            end
            fprintf('  pair TIPs used (%d pairs, min %.4f, max %.4f):\n', height(pr), ...
                    min(pr.TIP), max(pr.TIP));
            disp(pr);
            if bad
                warning('tip:NonNegativePairTIP', ['Correlation gives TIP >= 0 for at least one ' ...
                    'applied pair: %s / %s / M = %s REJECTED, predictions skipped.'], cc, md, Mdef);
                accRows(end+1, :) = {cc, md, Mdef, false, false, false, false, "rejected: TIP >= 0 for a pair"}; %#ok<AGROW>
                continue;
            end

            % --- binaries with the correlated TIP
            runner.clearWarm();
            errCor = nan(numel(binIdx), 1);  okCor = false(numel(binIdx), 1);
            shiftCor = nan(numel(binIdx), 1);  pctCor = nan(numel(binIdx), 1);
            for n = 1:numel(binIdx)
                i = binIdx(n);
                c = dp.getCase(i);
                g = find(cfg.groups == groupKey(i), 1);
                tC = tipCor(g);
                r  = runner.solve(c, md, [0, tC; tC, 0], 'Key', sprintf("cor|%s|%d|%s", Mdef, i, md));
                sh = shiftOf(r, cfg);
                errCor(n) = sh - dp.ShiftMid_psi(i);
                okCor(n)  = r.ok;
                shiftCor(n) = sh;
                pctCor(n) = pctErr(r.P, dp.ConfMid_psi(i), psi);
                binRows(end+1, :) = {cc, md, Mdef, c.Mixture, c.Rock, c.z(end), tipOpt(g), tC, ...
                    dp.ShiftMid_psi(i), dp.ConfHalfWidth_psi(i), sh, errCor(n), ...
                    abs(errCor(n)) <= dp.ConfHalfWidth_psi(i), pctCor(n), r.ok, r.route, r.reason}; %#ok<AGROW>
            end
            rmseCor = sqrt(mean(errCor(okCor).^2));

            % --- ternaries, blind
            terErr = [];  terIn = [];  terPct = [];
            for i = find(isTer).'
                c   = dp.getCase(i);
                eng = runner.engine(c.Mixture, c.Rock);
                K   = pairMatrix(eng, cfg, c, Mdef, tipOf);
                r   = runner.solve(c, md, K, 'Key', sprintf("ter|%s|%d|%s", Mdef, i, md));
                sh  = shiftOf(r, cfg);
                err = sh - dp.ShiftMid_psi(i);
                inB = abs(err) <= dp.ConfHalfWidth_psi(i);
                terErr(end+1) = err;  terIn(end+1) = inB && r.ok; %#ok<AGROW>
                terPct(end+1) = pctErr(r.P, dp.ConfMid_psi(i), psi); %#ok<AGROW>
                terRows(end+1, :) = {cc, md, Mdef, c.Mixture, c.Rock, dp.ShiftMid_psi(i), ...
                    dp.ConfHalfWidth_psi(i), sh, err, inB, terPct(end), r.ok, ...
                    r.route, r.reason}; %#ok<AGROW>
            end

            % --- gap and acceptance
            gapRows(end+1, :) = {cc, md, Mdef, role, rmseOpt, rmseCor, rmseCor - rmseOpt, ...
                nnz(okOpt), nnz(okCor), numel(binIdx), sqrt(mean(terErr.^2, 'omitnan'))}; %#ok<AGROW>
            isC5  = dp.Mixture(binIdx) == "C1-nC5";
            isC10 = dp.Mixture(binIdx) == "C1-nC10";
            a1 = ~cfg.accept.noNegativeC1nC5 || all(shiftCor(isC5) > 0);
            a2 = max(abs(errCor(isC10)), [], 'omitnan') <= cfg.accept.maxErr_C1nC10_psi;
            a3 = ~cfg.accept.ternaryInBand || all(terIn);
            a4 = max([abs(pctCor); abs(terPct(:))], [], 'omitnan') <= cfg.accept.maxPdewErr_pct;
            accRows(end+1, :) = {cc, md, Mdef, a1, a2, a3, a4, ...
                string(verdict(a1 && a2 && a3 && a4))}; %#ok<AGROW>
            fprintf(['  binary shift RMSE: optimised %.1f -> correlated %.1f psi (gap %+.1f) | ' ...
                     'ternary RMSE %.1f psi\n  criteria: C1-nC5 sign %d | C1-nC10 %d | ternaries %d | ' ...
                     '1%% Pdew %d -> %s\n'], rmseOpt, rmseCor, rmseCor - rmseOpt, ...
                    sqrt(mean(terErr.^2, 'omitnan')), a1, a2, a3, a4, verdict(a1 && a2 && a3 && a4));
        end
    end
end

%% 4. Tables and output --------------------------------------------------------
fitTable = cell2table(fitRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', 'Role', 'Form', ...
    'alpha', 'se_alpha', 'beta', 'se_beta', 'gamma', 'se_gamma', 'RMSE_TIP', 'RMSE_logfit', ...
    'dof', 'Leverage_B1', 'Note'});
levTable = cell2table(levRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', 'Role', 'Group', ...
    'sigma_over_rp', 'M', 'Leverage', 'TIP_optimised', 'TIP_correlated', 'Delta_TIP', ...
    'gamma_dropOne', 'dgamma'});
binTable = cell2table(binRows, 'VariableNames', {'Combination', 'Mode', 'TIPsource', 'Mixture', ...
    'Rock', 'y_heavy', 'TIP_optimised', 'TIP_used', 'Shift_exp_psi', 'HalfWidth_psi', ...
    'Shift_pred_psi', 'Error_psi', 'InBand', 'Pdew_err_pct', 'Accepted', 'Route', 'Reason'});
terTable = cell2table(terRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', 'Mixture', 'Rock', ...
    'Shift_exp_psi', 'HalfWidth_psi', 'Shift_pred_psi', 'Error_psi', 'InBand', 'Pdew_err_pct', ...
    'Accepted', 'Route', 'Reason'});
gapTable = cell2table(gapRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', 'Role', ...
    'RMSE_shift_optimised_psi', 'RMSE_shift_correlated_psi', 'Gap_psi', 'nAccepted_opt', ...
    'nAccepted_cor', 'nBinaries', 'RMSE_ternary_psi'});
pairTableAll = cell2table(pairRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', 'Context', ...
    'Rock', 'Pair', 'M', 'sigma_over_rp', 'TIP'});
accTable = cell2table(accRows, 'VariableNames', {'Combination', 'Mode', 'Mdef', ...
    'C1nC5_sign', 'C1nC10_error', 'Ternaries_in_band', 'Pdew_within_1pct', 'Verdict'});
for nm = ["fitTable", "levTable", "binTable", "terTable", "gapTable", "pairTableAll", "accTable"]
    t = eval(nm);
    isTxt = varfun(@(v) iscell(v) && all(cellfun(@(e) isstring(e) || ischar(e), v)), t, ...
                   'OutputFormat', 'uniform');
    if any(isTxt), t = convertvars(t, t.Properties.VariableNames(isTxt), 'string'); end
    eval(nm + " = t;");
end

fprintf('\n=== Correlation fits (%s, 3 parameters, 4 groups, dof = 1) ===\n', cfg.tipForm);
disp(fitTable);
fprintf('=== Leverage, correlated vs optimised TIP, drop-one gamma ===\n'); disp(levTable);
fprintf('=== Binary predictions ===\n');                                    disp(binTable);
fprintf('=== Ternary predictions (blind) ===\n');                           disp(terTable);
fprintf('=== Calibration-to-prediction gap ===\n');                         disp(gapTable);
fprintf('=== Pair TIPs used in every prediction ===\n');                    disp(pairTableAll);
fprintf('=== Acceptance against the pre-declared criteria ===\n');          disp(accTable);

save(fullfile(cfg.resultsDir, cfg.outFile + ".mat"), 'cfg', 'regTable', 'fitTable', 'levTable', ...
     'binTable', 'terTable', 'gapTable', 'pairTableAll', 'accTable');
xls = fullfile(cfg.resultsDir, cfg.outFile + ".xlsx");
if isfile(xls), delete(xls); end
writetable(regTable,     xls, 'Sheet', 'Regression_gate');
writetable(fitTable,     xls, 'Sheet', 'Correlation');
writetable(levTable,     xls, 'Sheet', 'Leverage');
writetable(binTable,     xls, 'Sheet', 'Binaries');
writetable(terTable,     xls, 'Sheet', 'Ternaries');
writetable(gapTable,     xls, 'Sheet', 'Gap');
writetable(pairTableAll, xls, 'Sheet', 'Pair_TIPs');
writetable(accTable,     xls, 'Sheet', 'Acceptance');
fprintf('\nWritten: %s\n', xls);

%% ============================== Local functions ==============================
function f = makeTipFun(par, form)
    if form == "loglinear"
        f = @(M, sR) -exp(par(1) + par(2) ./ M + par(3) * sR);
    else
        f = @(M, sR) par(1) + par(2) ./ M + par(3) * sR;
    end
end

function sh = shiftOf(r, cfg)
    if cfg.shiftBasis == "expBulk"
        sh = r.shiftExpBulk_psi;      % P_model - P_bulk,experimental (V31)
    else
        sh = r.shift_psi;             % P_model - P_bulk,model
    end
end

function e = pctErr(P_Pa, Pexp_psi, psi)
    e = 100 * (P_Pa / (Pexp_psi * psi) - 1);
end

function s = verdict(pass)
    if pass, s = "PASS"; else, s = "FAIL"; end
end

function K = pairMatrix(eng, cfg, c, Mdef, tipOf)
    % TIP matrix for a multicomponent case: every pair from the correlation.
    sig = eng.eos.Fluid.LJ_Size(:);
    MW  = eng.eos.Fluid.MW(:);
    rp  = cfg.poreRadius(char(c.Rock));
    nc  = numel(c.z);
    if Mdef == "heavy"
        M = max(MW(:), MW(:).');
    else
        M = sqrt(MW(:) * MW(:).');
    end
    sR = 0.5 * (sig(:) + sig(:).') / rp;
    K = tipOf(M, sR);
    K(1:nc+1:end) = 0;
    K = (K + K.') / 2;
end

function [T, bad] = pairTable(runner, dp, cfg, binIdx, isTer, Mdef, tipOf)
    % Every pair the correlation is applied to, with its TIP. bad = any TIP >= 0.
    rows = {};
    seen = strings(0, 1);
    add = @(ctx, rock, pair, M, sR) {ctx, rock, pair, M, sR, tipOf(M, sR)};
    for i = [binIdx; find(isTer)].'
        c   = dp.getCase(i);
        eng = runner.engine(c.Mixture, c.Rock);
        sig = eng.eos.Fluid.LJ_Size(:);
        MW  = eng.eos.Fluid.MW(:);
        rp  = cfg.poreRadius(char(c.Rock));
        nm  = eng.eos.Fluid.ComponentNames;
        ctx = "binary";
        if numel(c.z) > 2, ctx = "ternary"; end
        for a = 1:numel(c.z) - 1
            for b = a + 1:numel(c.z)
                key = ctx + "|" + c.Rock + "|" + nm(a) + "-" + nm(b);
                if any(seen == key), continue; end
                seen(end+1) = key; %#ok<AGROW>
                if Mdef == "heavy", M = max(MW(a), MW(b)); else, M = sqrt(MW(a) * MW(b)); end
                sR = 0.5 * (sig(a) + sig(b)) / rp;
                rows(end+1, :) = add(ctx, c.Rock, nm(a) + "-" + nm(b), M, sR); %#ok<AGROW>
            end
        end
    end
    T = cell2table(rows, 'VariableNames', {'Context', 'Rock', 'Pair', 'M', 'sigma_over_rp', 'TIP'});
    T = convertvars(T, {'Context', 'Rock', 'Pair'}, 'string');
    bad = any(T.TIP >= 0);
end