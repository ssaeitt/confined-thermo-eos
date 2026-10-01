%% TC_PERFORMANCE  Critical-shift performance of each locked parameter set
%
% For every (k, p, lambda) combination this reports the confined critical
% temperature shift under BOTH measures:
%   "exact"      the confined PR critical point (ConfinedEOS.pureCriticalPoint),
%                the same EOS that computes the dew points;
%   "linearized" the V31 expression dTc = c_i(Tc,conf,exp) / a_i(Tc), which sets
%                alpha = 1 and therefore overstates dTc by roughly (1 + m_i).
%
% Purpose: state joint-linearized's EXACT-EOS critical-shift performance
% explicitly. In a joint fit dominated by the mixture block, the linearized
% relation acts as a REGULARISING CONSTRAINT on (k, p, lambda) - it keeps the
% confinement term on a physically sensible scale - and must not be presented
% as the calibration of the model. The exact-EOS numbers below are what the
% model actually delivers on the 59 literature points.
%
% Calibration wall: quartz only (the literature dTc data are silica pores).

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
            error('tc:PackageNotFound', 'Folder "+%s" not found under %s.', pk, scriptDir);
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
cfg.wallRock    = "EF2";           % row used only to read mineral_E
cfg.wallMineral = "Quartz";        % calibration wall: quartz only
cfg.outFile     = "tc_performance";

cfg.combos = table( ...
    ["sequential-exact"; "joint-exact"; "joint-linearized"], ...
    [151.385346647943;  98.0578847776961;  112.42634826677 ], ...
    [1.30682834548441;  1.69745310322519;  1.3558179290678 ], ...
    [4.51381302592871;  10.2425163919516;  2.0916235090532 ], ...
    'VariableNames', {'Combination', 'k', 'p', 'lambda'});

if ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end

%% 1. Data and calibration engine ---------------------------------------------
fluid   = entities.FluidProperties.loadFromWorkbook(cfg.mixtureFile);
rockRef = entities.RockProperties.loadFromWorkbook(cfg.mixtureFile, cfg.wallRock);
q       = rockRef.mineral(cfg.wallMineral);
calRock = entities.RockProperties("CalWall-" + q.Name, 0, q.Name, q.Epsilon_K, q.GrainDensity, 1);
eos     = thermo.ConfinedEOS(fluid, calRock);
fprintf('Calibration wall: %s only, E = %.4f K (dTc data are silica pores).\n', ...
        cfg.wallMineral, q.Epsilon_K);

tcAll = entities.TcShiftDataset.loadFromWorkbook(cfg.tcShiftFile, 'Sheet', cfg.tcSheet);
tc    = tcAll.select(ismember(tcAll.Component, fluid.ComponentNames));
[idx, r] = tc.modelInputs(fluid);
fprintf('dTc points: %d (%s)\n\n', tc.NumPoints, strjoin(unique(tc.Component, 'stable'), ", "));

%% 2. Both measures for every combination -------------------------------------
sumRows = {};  compRows = {};  pointTables = struct();
for ci = 1:height(cfg.combos)
    comb = cfg.combos.Combination(ci);
    eos.k = cfg.combos.k(ci);  eos.pT_wall = cfg.combos.p(ci);  eos.lambda = cfg.combos.lambda(ci);

    dLin = predictTc(eos, tc, idx, r, "linearized");
    dEx  = predictTc(eos, tc, idx, r, "exact");
    sLin = metrics(dLin, tc.dTc);
    sEx  = metrics(dEx,  tc.dTc);

    sumRows(end+1, :) = {comb, cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci), ...
        sEx.rmse, sEx.aard, sEx.r2, sEx.bias, sLin.rmse, sLin.aard, sLin.r2, sLin.bias, ...
        sEx.n, tc.NumPoints}; %#ok<AGROW>

    for comp = unique(tc.Component, 'stable').'
        m = tc.Component == comp;
        cEx  = metrics(dEx(m),  tc.dTc(m));
        cLin = metrics(dLin(m), tc.dTc(m));
        compRows(end+1, :) = {comb, comp, nnz(m), cEx.rmse, cEx.aard, cEx.bias, ...
            cLin.rmse, cLin.aard, cLin.bias, mean(dEx(m) ./ dLin(m), 'omitnan')}; %#ok<AGROW>
    end

    pointTables.(matlab.lang.makeValidName(comb)) = [tc.toTable(), ...
        table(dEx, dLin, dEx - tc.dTc, dLin - tc.dTc, ...
        'VariableNames', {'dTc_exact', 'dTc_linearized', 'Res_exact', 'Res_linearized'})];

    fprintf(['%-18s k = %8.3f, p = %.4f, lambda = %7.4f\n' ...
             '   exact EOS  : RMSE %.4f | AARD %5.1f %% | R2 %+.3f | bias %+.4f\n' ...
             '   linearized : RMSE %.4f | AARD %5.1f %% | R2 %+.3f | bias %+.4f\n'], ...
            comb, cfg.combos.k(ci), cfg.combos.p(ci), cfg.combos.lambda(ci), ...
            sEx.rmse, sEx.aard, sEx.r2, sEx.bias, sLin.rmse, sLin.aard, sLin.r2, sLin.bias);
end

summaryTable = cell2table(sumRows, 'VariableNames', {'Combination', 'k', 'p', 'lambda', ...
    'RMSE_exact', 'AARD_exact_pct', 'R2_exact', 'Bias_exact', ...
    'RMSE_lin', 'AARD_lin_pct', 'R2_lin', 'Bias_lin', 'nUsed', 'nPoints'});
compTable = cell2table(compRows, 'VariableNames', {'Combination', 'Component', 'n', ...
    'RMSE_exact', 'AARD_exact_pct', 'Bias_exact', 'RMSE_lin', 'AARD_lin_pct', 'Bias_lin', ...
    'ratio_exact_over_lin'});
summaryTable = convertvars(summaryTable, {'Combination'}, 'string');
compTable = convertvars(compTable, {'Combination', 'Component'}, 'string');

fprintf('\n=== Critical-shift performance ===\n');  disp(summaryTable);
fprintf('=== Per component ===\n');                 disp(compTable);

%% 3. Statement for the paper ---------------------------------------------------
jl = summaryTable(summaryTable.Combination == "joint-linearized", :);
se = summaryTable(summaryTable.Combination == "sequential-exact", :);
statement = sprintf([ ...
    'joint-linearized reproduces the 59 literature critical-shift points with an EXACT-EOS ' ...
    'RMSE of %.4f (AARD %.1f %%, R2 %+.3f), against %.4f (AARD %.1f %%, R2 %+.3f) under the ' ...
    'linearized expression used inside its objective, and against %.4f (AARD %.1f %%) for the ' ...
    'sequential-exact parameters fitted on the exact criterion. In the joint fit the linearized ' ...
    'relation therefore acts as a regularising constraint on (k, p, lambda) in a ' ...
    'mixture-dominated objective; it is not the calibration of the model, and the exact-EOS ' ...
    'figures above are what the model delivers on critical shifts.'], ...
    jl.RMSE_exact, jl.AARD_exact_pct, jl.R2_exact, jl.RMSE_lin, jl.AARD_lin_pct, jl.R2_lin, ...
    se.RMSE_exact, se.AARD_exact_pct);
fprintf('\n=== Statement ===\n%s\n', statement);

save(fullfile(cfg.resultsDir, cfg.outFile + ".mat"), 'cfg', 'summaryTable', 'compTable', ...
     'pointTables', 'statement');
xls = fullfile(cfg.resultsDir, cfg.outFile + ".xlsx");
if isfile(xls), delete(xls); end
writetable(summaryTable, xls, 'Sheet', 'Summary');
writetable(compTable,    xls, 'Sheet', 'PerComponent');
writetable(table(string(statement), 'VariableNames', {'Statement'}), xls, 'Sheet', 'Statement');
fn = fieldnames(pointTables);
for i = 1:numel(fn)
    writetable(pointTables.(fn{i}), xls, 'Sheet', fn{i});
end
fprintf('\nWritten: %s\n', xls);

%% ============================== Local functions ==============================
function pred = predictTc(eos, tc, idx, r, model)
    pred = nan(tc.NumPoints, 1);
    for n = 1:tc.NumPoints
        i  = idx(n);
        Tc = eos.Fluid.Tc(i);
        try
            switch model
                case "linearized"
                    [ac, ~, ~] = eos.pureParameters(Tc, Inf);
                    [~, ~, c]  = eos.pureParameters(Tc * (1 - tc.dTc(n)), r(n));
                    pred(n) = c(i) / ac(i);
                case "exact"
                    pred(n) = 1 - eos.pureCriticalPoint(i, r(n)) / Tc;
            end
        catch ME
            if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
        end
    end
end

function s = metrics(pred, obs)
    ok = isfinite(pred) & isfinite(obs);
    e  = pred(ok) - obs(ok);
    s = struct('n', nnz(ok), 'rmse', sqrt(mean(e.^2)), ...
               'aard', 100 * mean(abs(e ./ obs(ok))), 'bias', mean(e), ...
               'r2', 1 - sum(e.^2) / sum((obs(ok) - mean(obs(ok))).^2));
end