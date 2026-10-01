%% PARETO_WEIGHT_SWEEP  Joint refit on the EXACT critical-shift criterion with
%  normalised residual blocks and a swept T_c weight.
%
%   J(w) = w * RMSE_dTc / s_dTc  +  RMSE_shift / s_shift
%
% Both blocks are divided by the RMS of their own experimental values, so the
% weight w is the only lever and w = 1 means "equal relative emphasis". This
% replaces the ad-hoc 0.05 / 100 scaling of the V31 objective, under which the
% exact-Tc joint fit is dominated by the linearized one.
%
% Decision vector: x = [k, p, lambda, TIP(1..G)] for ONE mode (cfg.mode).
% dTc uses the exact confined critical point; the mixture block uses the binary
% confined dew points through solvers.DewPointRunner, with the same acceptance
% guards as everywhere else. Rejected points are penalised in the objective but
% never enter the reported statistics.
%
% Output: the Pareto front of exact dTc RMSE against mixture shift RMSE, with
% the joint-linearized and sequential-exact reference points marked, plus a
% dominance test: does a weighted exact fit reach joint-linearized's mixture
% performance with a consistent critical-shift model?

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
            error('pareto:PackageNotFound', 'Folder "+%s" not found under %s.', pk, scriptDir);
        end
        addpath(fileparts(cand(1)));
    end
end

%% 0. Configuration ---------------------------------------------------------
cfg.configDir   = fullfile(scriptDir, "config");
cfg.resultsDir  = fullfile(scriptDir, "results");
cfg.mixtureFile = fullfile(cfg.configDir, "MixtureData.xlsx");
cfg.tcShiftFile = fullfile(cfg.configDir, "Extended_Dataset - Original.xlsx");
cfg.tcSheet     = "All Data";
cfg.outFile     = "pareto_weight_sweep";

cfg.T          = 293.15;
cfg.poreRadius = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
cfg.mode       = "FWIOnly";            % one mode per sweep (Combined differs by < 1 psi here)
cfg.wallRock   = "EF2";
cfg.wallMineral = "Quartz";            % calibration wall for dTc only

cfg.weights    = [0.1 0.3 0.5 1 2 3 10 30 100];   % T_c weight w (30, 100 show the
                                                  % sequential-exact limit)
cfg.variants   = ["exact"m, "linearized"];  % dTc measure INSIDE the objective;
                                           % both fronts are reported on the EXACT metric
cfg.groups     = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"];
cfg.x0         = [151.385346647943, 1.30682834548441, 4.51381302592871, ...
                  -3.48386287716462, -1.46205964046861, -1.12088738883485, -0.905508761783693];
cfg.lb         = [0.1, 0.1,  0, -100, -100, -100, -100];
cfg.ub         = [300, 2.5, 20,   50,   50,   50,   50];
cfg.maxIter    = 150;
cfg.fdStep     = 1e-6;
cfg.maxRestarts = 4;                   % restarts until ExitFlag 1 and no parameter at a bound
cfg.boundTolRel = 1e-3;                % "at a bound" if within this fraction of the range
cfg.autoExpandBounds = true;           % double a binding bound (within hard limits) and retry
cfg.hardLB     = [0.01, 0.05,  0, -1e4, -1e4, -1e4, -1e4];
cfg.hardUB     = [5e3,  4.0,  200,  200,  200,  200,  200];
cfg.carryWeights = [1 3];              % joint-exact points carried through the pipeline
cfg.carryLinearizedBest = true;        % plus the best linearized-objective point, if competitive
cfg.selectedFile = "pareto_selected.csv";  % read by TIP_correlation.m
cfg.modes_out  = ["FWIOnly", "Combined"];  % modes written for the carried points
                                           % (the sweep itself runs cfg.mode)
cfg.penalty    = 1e3;                  % psi, per rejected dew point (objective only)
cfg.warmStart  = true;                 % start each w from the previous solution

% Reference points to mark on the front (from the locked runs)
cfg.refs = table(["joint-linearized"; "sequential-exact"; "joint-exact"], ...
    [112.42634826677; 151.385346647943; 98.05788478], ...
    [1.3558179290678; 1.30682834548441; 1.69745310322519], ...
    [2.0916235090532; 4.51381302592871; 10.2425163919516], ...
    [-4.9229696302796,  -1.96920366521686, -1.44628069292089, -0.754862049795788; ...
     -3.48386287716462, -1.46205964046861, -1.12088738883485, -0.905508761783693; ...
     -13.4437287153297, -5.90838752675671, -4.59334889956389, -2.38867380361359], ...
    'VariableNames', {'Combination', 'k', 'p', 'lambda', 'TIP'});

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
psi = entities.DewPointDataset.psi2Pa;

%% 1. Data, engines, calibration wall ------------------------------------------
fluidAll = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile);
rockRef  = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallRock);
q        = rockRef.mineral(cfg.wallMineral);
calRock  = entities.RockProperties("CalWall-" + q.Name, 0, q.Name, q.Epsilon_K, q.GrainDensity, 1);
eosCal   = thermo.ConfinedEOS(fluidAll, calRock);

tcAll = entities.TcShiftDataset.loadFromWorkbook(cfg.tcShiftFile, 'Sheet', cfg.tcSheet);
tc    = tcAll.select(ismember(tcAll.Component, fluidAll.ComponentNames));
[tcIdx, tcR] = tc.modelInputs(fluidAll);

dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
groupKey = dp.Mixture + "|" + dp.Rock;
isBin  = dp.HasConfined & dp.NumComponents == 2 & ismember(groupKey, cfg.groups);
binIdx = find(isBin);
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius);
for i = binIdx.'
    runner.engine(dp.Mixture(i), dp.Rock(i));
    [~, ~, okb] = runner.bulk(dp.getCase(i));
    if ~okb, error('pareto:BulkFailed', 'Bulk dew point failed for case %d.', i); end
end
caseGroup = arrayfun(@(i) find(cfg.groups == groupKey(i), 1), binIdx);

% Block scales: RMS of the experimental values, so both blocks are ~O(1)
sTc    = sqrt(mean(tc.dTc.^2));
sShift = sqrt(mean(dp.ShiftMid_psi(binIdx).^2));
fprintf('Block scales: s_dTc = %.4f (%d points), s_shift = %.1f psi (%d points)\n', ...
        sTc, tc.NumPoints, sShift, numel(binIdx));
fprintf('Objective: J(w) = w * RMSE_dTc / s_dTc + RMSE_shift / s_shift, mode %s\n\n', cfg.mode);

ctx = struct('cfg', cfg, 'eosCal', eosCal, 'tc', tc, 'tcIdx', tcIdx, 'tcR', tcR, ...
             'dp', dp, 'binIdx', binIdx, 'caseGroup', caseGroup, 'runner', runner, ...
             'sTc', sTc, 'sShift', sShift);

%% 2. Reference points ----------------------------------------------------------
refRows = {};
for i = 1:height(cfg.refs)
    x = [cfg.refs.k(i), cfg.refs.p(i), cfg.refs.lambda(i), cfg.refs.TIP(i, :)];
    m = evaluate(x, ctx, "exact");
    refRows(end+1, :) = {cfg.refs.Combination(i), NaN, "reference", x(1), x(2), x(3), x(4), x(5), ...
        x(6), x(7), m.rmseTc, m.rmseTcLin, m.aardTc, m.rmseShift, m.nRejected, m.J1}; %#ok<AGROW>
    fprintf('reference %-18s exact dTc RMSE %.4f | shift RMSE %6.1f psi | rejected %d\n', ...
            cfg.refs.Combination(i), m.rmseTc, m.rmseShift, m.nRejected);
end

%% 3. Weight sweep, both objective variants -------------------------------------
sweepRows = {};
opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxIterations', cfg.maxIter, 'FiniteDifferenceStepSize', cfg.fdStep, ...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10);
for variant = cfg.variants
    fprintf('\n########## objective dTc metric: %s ##########\n', variant);
    x0 = cfg.x0;
    lb = cfg.lb;  ub = cfg.ub;
    for w = cfg.weights
        tW = tic;
        [x, flag, lb, ub, nRestart, atB] = solveConverged(@(v) objective(v, ctx, w, variant), ...
                                                x0, lb, ub, cfg, opts);
        m = evaluate(x, ctx, variant);
        okRun = (flag == 1) && ~any(atB);
        sweepRows(end+1, :) = {variant, w, x(1), x(2), x(3), x(4), x(5), x(6), x(7), ...
            m.rmseTc, m.rmseTcLin, m.aardTc, m.rmseShift, m.nRejected, m.J1, flag, ...
            strjoin(paramNames(atB), ","), nRestart, okRun, toc(tW) / 60}; %#ok<AGROW>
        fprintf(['%-10s w = %-6g exact dTc RMSE %.4f | shift RMSE %6.1f psi | rejected %d | ' ...
                 'flag %d | atBound [%s] | restarts %d | %.1f min\n'], variant, w, m.rmseTc, ...
                m.rmseShift, m.nRejected, flag, strjoin(paramNames(atB), ","), nRestart, toc(tW) / 60);
        if ~okRun
            warning('pareto:NotConverged', ['%s / w = %g did not meet the convergence standard ' ...
                '(ExitFlag %d, at bound: %s). Point flagged and excluded from the front.'], ...
                variant, w, flag, strjoin(paramNames(atB), ","));
        end
        if cfg.warmStart, x0 = x; end
    end
end

frontTable = cell2table(sweepRows, 'VariableNames', {'Variant', 'w', 'k', 'p', 'lambda', ...
    'TIP_C1nC5_EF2', 'TIP_C1nC8_EF2', 'TIP_C1nC10_EF2', 'TIP_C1nC8_B1', ...
    'RMSE_dTc_exact', 'RMSE_dTc_lin', 'AARD_dTc_pct', 'RMSE_shift_psi', 'nRejected', 'J_w1', ...
    'ExitFlag', 'AtBound', 'Restarts', 'Converged', 'Minutes'});
frontTable = convertvars(frontTable, {'Variant', 'AtBound'}, 'string');

refTable = cell2table(refRows, 'VariableNames', {'Combination', 'w', 'Variant', 'k', 'p', 'lambda', ...
    'TIP_C1nC5_EF2', 'TIP_C1nC8_EF2', 'TIP_C1nC10_EF2', 'TIP_C1nC8_B1', ...
    'RMSE_dTc_exact', 'RMSE_dTc_lin', 'AARD_dTc_pct', 'RMSE_shift_psi', 'nRejected', 'J_w1'});
refTable = convertvars(refTable, {'Combination', 'Variant'}, 'string');

% Pareto filter within each variant, converged points only
frontTable.OnParetoFront = false(height(frontTable), 1);
for v = cfg.variants
    sel = find(frontTable.Variant == v & frontTable.Converged);
    pts = [frontTable.RMSE_dTc_exact(sel), frontTable.RMSE_shift_psi(sel)];
    for i = 1:numel(sel)
        frontTable.OnParetoFront(sel(i)) = ~any(all(pts <= pts(i, :), 2) & any(pts < pts(i, :), 2));
    end
end
jl = refTable(refTable.Combination == "joint-linearized", :);
se = refTable(refTable.Combination == "sequential-exact", :);
frontTable.BeatsJointLinearized = frontTable.Converged & ...
    frontTable.RMSE_shift_psi <= jl.RMSE_shift_psi & frontTable.RMSE_dTc_exact <= jl.RMSE_dTc_exact;

fprintf('\n=== Reference points ===\n');  disp(refTable);
fprintf('=== Weight sweep (both variants, exact metric) ===\n');  disp(frontTable);
nBad = nnz(~frontTable.Converged);
if nBad > 0
    fprintf(2, '%d of %d sweep points failed the convergence standard; see AtBound / ExitFlag.\n', ...
            nBad, height(frontTable));
end
if any(frontTable.BeatsJointLinearized)
    b = frontTable(find(frontTable.BeatsJointLinearized, 1), :);
    fprintf(['\nVerdict: %s at w = %g dominates the original joint-linearized point ' ...
             '(dTc RMSE %.4f <= %.4f, shift RMSE %.1f <= %.1f psi).\n'], b.Variant, b.w, ...
            b.RMSE_dTc_exact, jl.RMSE_dTc_exact, b.RMSE_shift_psi, jl.RMSE_shift_psi);
else
    fprintf(['\nVerdict: no swept point dominates the original joint-linearized point; report ' ...
             'the trade-off curve and state the weight chosen.\n']);
end

%% 4. One figure, both fronts, endpoints marked ---------------------------------
figure('Color', 'w', 'Name', 'Pareto fronts');
hold on; grid on; box on;
sty = ["o-", "s--"];
for vi = 1:numel(cfg.variants)
    sel = frontTable.Variant == cfg.variants(vi) & frontTable.Converged;
    [~, o] = sort(frontTable.w(sel));
    T = frontTable(sel, :);  T = T(o, :);
    plot(T.RMSE_dTc_exact, T.RMSE_shift_psi, sty(vi), 'LineWidth', 1.5, ...
         'DisplayName', "joint-" + cfg.variants(vi) + " objective");
    text(T.RMSE_dTc_exact, T.RMSE_shift_psi, compose("  w=%g", T.w), 'FontSize', 8);
    bad = frontTable.Variant == cfg.variants(vi) & ~frontTable.Converged;
    if any(bad)
        plot(frontTable.RMSE_dTc_exact(bad), frontTable.RMSE_shift_psi(bad), 'rx', ...
             'MarkerSize', 9, 'LineWidth', 1.2, 'DisplayName', cfg.variants(vi) + " (not converged)");
    end
end
plot(se.RMSE_dTc_exact, se.RMSE_shift_psi, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', [0.9 0.7 0.1], ...
     'DisplayName', 'sequential-exact (large-w endpoint)');
plot(jl.RMSE_dTc_exact, jl.RMSE_shift_psi, 'ks', 'MarkerSize', 11, 'MarkerFaceColor', [0.4 0.4 0.4], ...
     'DisplayName', 'original joint-linearized');
je = refTable(refTable.Combination == "joint-exact", :);
plot(je.RMSE_dTc_exact, je.RMSE_shift_psi, 'kd', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
     'DisplayName', 'original joint-exact');
xlabel('exact-EOS \Delta T_c RMSE');  ylabel('binary shift RMSE (psi)');
title(sprintf('Pareto fronts, mode %s (both objectives evaluated on the exact metric)', cfg.mode));
legend('Location', 'best');
saveas(gcf, fullfile(cfg.resultsDir, cfg.outFile + "_fronts.png"));

%% 5. Points carried into the pipeline -------------------------------------------
carry = {};
for w = cfg.carryWeights
    r = frontTable(frontTable.Variant == "exact" & frontTable.w == w & frontTable.Converged, :);
    if isempty(r)
        warning('pareto:CarryMissing', 'No converged exact point at w = %g to carry.', w);
        continue;
    end
    carry(end+1, :) = {sprintf("pareto-exact-w%g", w), r}; %#ok<AGROW>
end
if cfg.carryLinearizedBest
    sel = frontTable(frontTable.Variant == "linearized" & frontTable.Converged, :);
    if ~isempty(sel)
        [~, j] = min(sel.RMSE_shift_psi ./ sShift + sel.RMSE_dTc_exact ./ sTc);
        carry(end+1, :) = {sprintf("pareto-lin-w%g", sel.w(j)), sel(j, :)}; %#ok<AGROW>
    end
end

selRows = cell(0,9);
for i = 1:size(carry, 1)
    r = carry{i, 2};
    for md = cfg.modes_out
        selRows(end+1, :) = {carry{i, 1}, md, r.k, r.p, r.lambda, r.TIP_C1nC5_EF2, ...
            r.TIP_C1nC8_EF2, r.TIP_C1nC10_EF2, r.TIP_C1nC8_B1}; %#ok<AGROW>
    end
end
selectedTable = cell2table(selRows, 'VariableNames', {'Combination', 'Mode', 'k', 'p', 'lambda', ...
    'TIP1', 'TIP2', 'TIP3', 'TIP4'});
selectedTable = convertvars(selectedTable, {'Combination', 'Mode'}, 'string');
selFile = fullfile(cfg.resultsDir, cfg.selectedFile);
writetable(selectedTable, selFile);
fprintf('\nCarried points written to %s:\n', selFile);
disp(selectedTable);
fprintf(['Run TIP_correlation.m next: it picks this file up automatically (cfg.extraCombosFile) ' ...
         'and runs the correlation, the binaries with correlated TIP, both ternaries and the ' ...
         '1%% check for these parameter sets.\n']);

save(fullfile(cfg.resultsDir, cfg.outFile + ".mat"), 'cfg', 'frontTable', 'refTable', 'selectedTable');
xls = fullfile(cfg.resultsDir, cfg.outFile + ".xlsx");
if isfile(xls), delete(xls); end
writetable(frontTable,    xls, 'Sheet', 'Front');
writetable(refTable,      xls, 'Sheet', 'References');
writetable(selectedTable, xls, 'Sheet', 'Carried');
fprintf('\nWritten: %s\n', xls);

%% ============================== Local functions ==============================
function nm = paramNames(mask)
    all = ["k", "p", "lambda", "TIP1", "TIP2", "TIP3", "TIP4"];
    nm = all(mask);
    if isempty(nm), nm = strings(1, 0); end
end

function [x, flag, lb, ub, nRestart, atB] = solveConverged(fun, x0, lb, ub, cfg, opts)
    % fmincon until ExitFlag 1 with no parameter at a bound. A binding bound is
    % doubled (within cfg.hardLB / cfg.hardUB) and the solve restarted, so the
    % reported optimum is interior rather than an artefact of the box.
    nRestart = 0;
    x = x0;
    while true
        [x, ~, flag] = fmincon(fun, x, [], [], [], [], lb, ub, [], opts);
        rng_ = ub - lb;
        atB = (x - lb) <= cfg.boundTolRel * rng_ | (ub - x) <= cfg.boundTolRel * rng_;
        if flag == 1 && ~any(atB), return; end
        if nRestart >= cfg.maxRestarts, return; end
        nRestart = nRestart + 1;
        if any(atB) && cfg.autoExpandBounds
            loB = (x - lb) <= cfg.boundTolRel * rng_;
            hiB = (ub - x) <= cfg.boundTolRel * rng_;
            lb(loB) = max(cfg.hardLB(loB), lb(loB) - abs(rng_(loB)));
            ub(hiB) = min(cfg.hardUB(hiB), ub(hiB) + abs(rng_(hiB)));
            x = min(max(x, lb), ub);
        else
            x = min(max(x .* (1 + 0.02 * randn(size(x))), lb), ub);   % nudge and retry
        end
    end
end
function J = objective(x, ctx, w, tcMetric)
    m = evaluate(x, ctx, tcMetric);
    J = w * m.rmseTcObj / ctx.sTc + m.rmseShiftPenalised / ctx.sShift;
end

function m = evaluate(x, ctx, tcMetric)
    % dTc block (both measures) and binary-shift block for x = [k p lambda TIPs].
    % tcMetric selects which dTc measure enters the OBJECTIVE; the reported
    % RMSE_dTc is always the exact one.
    if nargin < 3, tcMetric = "exact"; end
    ctx.eosCal.k = x(1);  ctx.eosCal.pT_wall = x(2);  ctx.eosCal.lambda = x(3);
    [predEx, predLin] = deal(nan(ctx.tc.NumPoints, 1));
    for n = 1:ctx.tc.NumPoints
        i  = ctx.tcIdx(n);
        Tc = ctx.eosCal.Fluid.Tc(i);
        try
            predEx(n) = 1 - ctx.eosCal.pureCriticalPoint(i, ctx.tcR(n)) / Tc;
        catch ME
            if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
        end
        try
            [ac, ~, ~] = ctx.eosCal.pureParameters(Tc, Inf);
            [~, ~, cc] = ctx.eosCal.pureParameters(Tc * (1 - ctx.tc.dTc(n)), ctx.tcR(n));
            predLin(n) = cc(i) / ac(i);
        catch ME
            if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
        end
    end
    eEx = predEx - ctx.tc.dTc;   eEx(~isfinite(eEx)) = 1;
    eLin = predLin - ctx.tc.dTc; eLin(~isfinite(eLin)) = 1;
    m.rmseTc    = sqrt(mean(eEx.^2));              % reported metric: exact, always
    m.rmseTcLin = sqrt(mean(eLin.^2));
    m.rmseTcObj = m.rmseTc;
    if tcMetric == "linearized", m.rmseTcObj = m.rmseTcLin; end
    ok = isfinite(predEx);
    m.aardTc = 100 * mean(abs(eEx(ok) ./ ctx.tc.dTc(ok)));

    ctx.runner.setFWI(x(1), x(2), x(3));
    tips = x(4:end);
    res = nan(numel(ctx.binIdx), 1);
    nRej = 0;
    for n = 1:numel(ctx.binIdx)
        i = ctx.binIdx(n);
        c = ctx.dp.getCase(i);
        t = tips(ctx.caseGroup(n));
        r = ctx.runner.solve(c, ctx.cfg.mode, [0, t; t, 0], ...
                             'Key', sprintf("sweep|%d|%s", i, ctx.cfg.mode));
        if r.ok
            res(n) = r.shift_psi - ctx.dp.ShiftMid_psi(i);
        else
            nRej = nRej + 1;
        end
    end
    m.nRejected = nRej;
    m.rmseShift = sqrt(mean(res.^2, 'omitnan'));                    % accepted only (statistics)
    m.rmseShiftPenalised = sqrt((sum(res.^2, 'omitnan') + nRej * ctx.cfg.penalty^2) / numel(res));
    m.J1 = m.rmseTc / ctx.sTc + m.rmseShift / ctx.sShift;           % w = 1 reference value
end