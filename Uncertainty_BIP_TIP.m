%% UNCERTAINTY_BIP_TIP  Measurement-based BIP windows and TIP corner sets
%
% Replaces the fixed +-1 % bulk-Pdew rule by a per-point band built from the
% laboratory accuracy figures, then re-derives the BIP windows and the TIP
% corner sets (rock-offset correlation) from those windows.
%
% Laboratory figures (experimental section, visual PVT cell):
%   temperature    5 % of the reading in degC (+-1 degC at 20 degC, +-2 at 40, +-3 at 60)
%   composition    5-10 % relative on the heavy component(s) as charged
%   dew detection  by volume: first liquid read at 0.1 cc in a 20 cc cell
%                  (0.5 % of V_sat) -> measured Pdew lies BELOW the true one
%   liquid volume  read to 0.1 cc -> MLDO +-0.5 points
%   pressure gauge high accuracy (not propagated); 1 psi/min depletion, equilibrium assumed
%
% 1. Per-point budget (solvers.UncertaintyBudget, bulk EOS, workbook BIPs):
%      UT  = max |Pdew(T +- dT) - Pdew|
%      UX  = RSS over heavy components of max |Pdew(x_h (1 +- r)) - Pdew|, r = 5 % and 10 %
%      lag = Pdew - P(liquid volume = 0.5 % V_sat)               (one-sided)
%    Band for the model bulk Pdew:   [Pm - U, Pm + lag + U],  U = sqrt(UX^2 + UT^2)
%    Weighted LSQ: target Pm + lag/2, sigma^2 = (UX^2 + UT^2)/3 + lag^2/12
%                  (bounds read as rectangular distributions)
%    MLDO budget (cfg.mldo, sheet Bulk_PVT, all fluids and temperatures): the same
%    composition and temperature perturbations on the CCE MLDO, plus the 0.5-point
%    reading resolution; model - measured compared with U = sqrt(UX^2 + UT^2 + Uread^2).
%    Composition check: relative change of the heavy fraction(s) of every adjusted
%    Mix_Props composition against the nominal one, compared with the 5 / 10 % bound.
% 2. BIP window per pair: BIPs for which EVERY 20 degC hydrocarbon row containing
%    the pair (binaries, ternary, late-core bulk) lies inside its band; other BIPs
%    at the workbook values. Also reported: weighted-LSQ BIP and the old +-1 % window.
% 3. TIP corner sets: each group TIP refitted (FWI-Only, k, p, lambda frozen at jL)
%    at the low / centre / high BIP of its own pair (early-campaign rows, as the presets); the 27 corners combine them and
%    are each fitted by the rock-offset correlation (thermo.RockOffsetTIP.custom).
%    Output: corner table, alpha/B ranges, and (cfg.predict) the confined-Pdew
%    prediction envelope over the 8 extreme corners plus the centre.
%
% Needs (local): config/<workbook of cfg.variant>.  Nothing here changes the
% workbook or the frozen presets.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -----------------------------------------------------------------
cfg.variant = "recal";         % "recal" (Recalibrated, jL2; current model) | "frozen" (Complete Mixtures, jL)
variants = struct( ...
    'frozen', struct('file', "MixtureData - Complete Mixtures.xlsx", 'preset', "jL"), ...
    'recal',  struct('file', "MixtureData - Recalibrated.xlsx",     'preset', "jL2"));
vsel = variants.(char(cfg.variant));
cfg.file   = fullfile(scriptDir, "config", vsel.file);
cfg.preset = vsel.preset;                          % source of k, p, lambda
cfg.outDir = fullfile(scriptDir, "results", "uncertainty_bip_tip_" + cfg.variant);
cfg.T      = 293.15;

% Laboratory figures
cfg.lab = struct('CompRel', [0.05, 0.10], 'TempRel', 0.05, 'CellVolume_cc', 20, 'DetectVolume_cc', 0.1, ...
                 'Repeat_psi', 0);
cfg.mldo = true;                                  % MLDO budget + composition-adjustment check
% Repeat_psi: optional repeatability bound added in quadrature to U (0 = the
% laboratory figures only). The two C1-nC8 [0.9, 0.1] bulk values (4110 and 4126 psi,
% same composition) differ by 16 psi, i.e. a bound of about 8 psi.

% Pairs, BIP search ranges and the rows used
cfg.pairs     = {["C1", "nC5"], ["C1", "nC8"], ["C1", "nC10"]};
cfg.bipRange  = struct('C1_nC5', [-0.06, 0.06], 'C1_nC8', [-0.02, 0.06], 'C1_nC10', [0.0, 0.08]);
cfg.excludeComponents = "CO2";                    % CO2 rows are outside this study

% Confinement
cfg.poreRadius   = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
cfg.lateCase     = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.lateRadius_m = 7.82711826617637e-9;           % back-calculated late-core radius

% TIP corner sets
cfg.cornerScenarios = [0.05, 0.10];               % composition scenarios for the corners
cfg.cornerCentre    = "workbook";                 % "workbook" | "wls" (centre level of each pair)
cfg.tipTarget       = "shift";                    % "shift" (agreed) | "absolute"
cfg.tipSearch       = struct('grid', linspace(-6, 2.5, 35), 'tolX', 1e-6);   % ln(-TIP)
cfg.descriptors     = ["inv_wbar", "inv_Tstar"];
cfg.predict         = true;                       % prediction envelope (extreme corners + centre)

if ~isfile(cfg.file), error('unc:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.file); end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;
xls = fullfile(cfg.outDir, "uncertainty_bip_tip.xlsx");
if isfile(xls), delete(xls); end
rock = entities.RockProperties.loadFromWorkbook(cfg.file, "EF2");
fwiRef = thermo.RockOffsetTIP.frozen(cfg.preset, cfg.file);
ub = solvers.UncertaintyBudget('CompRel', cfg.lab.CompRel, 'TempRel', cfg.lab.TempRel, ...
        'CellVolume_cc', cfg.lab.CellVolume_cc, 'DetectVolume_cc', cfg.lab.DetectVolume_cc);
nS = numel(cfg.lab.CompRel);

%% 1. Bulk rows (20 degC, hydrocarbons, unique mixture + z + bulk Pdew) ---------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.file, 'DefaultT', cfg.T);
rows = struct('Mixture', {}, 'z', {}, 'T', {}, 'Pm', {}, 'Label', {});
for i = 1:dp.NumCases
    c = dp.getCase(i);
    if abs(c.T - cfg.T) > 1e-6 || any(ismember(c.Components, cfg.excludeComponents)), continue; end
    dup = arrayfun(@(r) r.Mixture == c.Mixture && numel(r.z) == numel(c.z) && ...
                   max(abs(r.z - c.z)) < 1e-9 && abs(r.Pm - c.PdewBulk) < 1e-6, rows);
    if any(dup), continue; end
    lab = c.Mixture + " [" + strjoin(compose("%.4g", c.z.'), ",") + "] " + ...
          string(round(c.PdewBulk / psi)) + " psi";
    rows(end+1) = struct('Mixture', c.Mixture, 'z', c.z, 'T', c.T, 'Pm', c.PdewBulk, 'Label', lab); %#ok<SAGROW>
end
nR = numel(rows);
fprintf('%d bulk rows at %.2f degC\n', nR, cfg.T - 273.15);

%% 2. Per-point uncertainty budget -------------------------------------------------------
U = repmat(struct('P0', NaN, 'UT', NaN, 'UX', nan(1, nS), 'lag', NaN, 'oneSided', false(1, nS), 'note', ""), nR, 1);
budRows = {};
for n = 1:nR
    fl = entities.FluidProperties.loadFromWorkbook(cfg.file, rows(n).Mixture);
    u  = ub.evaluate(fl, rock, rows(n).z, rows(n).T, rows(n).Pm);
    U(n) = struct('P0', u.P0, 'UT', u.UT, 'UX', u.UX, 'lag', u.lag, 'oneSided', u.oneSided, 'note', u.note);
    rec = {rows(n).Label, rows(n).Pm / psi, u.P0 / psi, 100 * (u.P0 / rows(n).Pm - 1), u.dP_T(1) / psi, ...
           u.dP_T(2) / psi, u.UT / psi};
    for q = 1:nS
        Uq = hypot(u.UX(q), u.UT);
        rec = [rec, {u.UX(q) / psi, Uq / psi, 100 * Uq / rows(n).Pm, u.oneSided(q)}]; %#ok<AGROW>
    end
    budRows(end+1, :) = [rec, {u.lag / psi, 100 * u.lag / rows(n).Pm, u.note}]; %#ok<SAGROW>
    fprintf('  %-42s UT %5.1f | UX %s | lag %5.1f psi %s\n', rows(n).Label, u.UT / psi, ...
            strjoin(compose("%6.1f", u.UX / psi), " / "), u.lag / psi, u.note);
end
vn = {'Row', 'Pmeas_psi', 'Pmodel_psi', 'Bulk_err_pct', 'dP_Tm1_psi', 'dP_Tp1_psi', 'UT_psi'};
for q = 1:nS
    tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
    vn = [vn, {['UX' tag '_psi'], ['U' tag '_psi'], ['U' tag '_pct'], ['UX' tag '_oneSided']}]; %#ok<AGROW>
end
vn = [vn, {'lag_psi', 'lag_pct', 'Note'}];
Budget = cell2table(budRows, 'VariableNames', vn);
Budget = convertvars(Budget, {'Row', 'Note'}, 'string');
writetable(Budget, xls, 'Sheet', 'Budget');
disp(Budget);

%% 2b. MLDO budget (Bulk_PVT: all fluids and temperatures) ----------------------------------
if cfg.mldo
    bp = readtable(cfg.file, 'Sheet', 'Bulk_PVT', 'VariableNamingRule', 'preserve');
    mRows = {};
    for i = 1:height(bp)
        mix = string(bp.Mixture(i));
        if any(ismember(split(mix, "-"), cfg.excludeComponents)) || ~isfinite(bp.MLDO_pct(i)), continue; end
        z  = str2double(split(erase(string(bp.mole_frac(i)), ["[", "]", " "]), ","));
        fl = entities.FluidProperties.loadFromWorkbook(cfg.file, mix);
        m  = ub.evaluateMLDO(fl, rock, z, bp.T(i), bp.Bulk_Pdew(i) * psi);
        d  = m.M0 - bp.MLDO_pct(i);
        rec = {string(bp.Fluid_ID(i)), mix, bp.T(i) - 273.15, m.dT, bp.MLDO_pct(i), m.M0, d, m.UT, m.Uread};
        for q = 1:nS
            Uq = sqrt(m.UX(q)^2 + m.UT^2 + m.Uread^2);
            rec = [rec, {m.UX(q), Uq, abs(d) <= Uq}]; %#ok<AGROW>
        end
        mRows(end+1, :) = [rec, {m.note}]; %#ok<SAGROW>
        fprintf('  MLDO %-4s %4.0f C: meas %4.1f | model %5.2f (%+5.2f) | UT %4.2f | UX %s | read %3.1f points %s\n', ...
                string(bp.Fluid_ID(i)), bp.T(i) - 273.15, bp.MLDO_pct(i), m.M0, d, m.UT, ...
                strjoin(compose("%5.2f", m.UX), " / "), m.Uread, m.note);
    end
    vnM = {'Fluid_ID', 'Mixture', 'T_C', 'dT_C', 'MLDO_meas_pct', 'MLDO_model_pct', 'Diff_pts', 'UT_pts', 'Uread_pts'};
    for q = 1:nS
        tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
        vnM = [vnM, {['UX' tag '_pts'], ['U' tag '_pts'], ['Inside' tag]}]; %#ok<AGROW>
    end
    MLDOBudget = cell2table(mRows, 'VariableNames', [vnM, {'Note'}]);
    MLDOBudget = convertvars(MLDOBudget, {'Fluid_ID', 'Mixture', 'Note'}, 'string');
    writetable(MLDOBudget, xls, 'Sheet', 'MLDO_Budget');
    disp(MLDOBudget);

    % Composition check: adjusted Mix_Props compositions vs the nominal (as charged) ones
    mp = readtable(cfg.file, 'Sheet', 'Mix_Props', 'VariableNamingRule', 'preserve');
    if all(ismember(["Nominal_mole_frac", "Fluid_ID"], string(mp.Properties.VariableNames)))
        aRows = {};  seen = strings(0, 1);
        for i = 1:height(mp)
            fid = string(mp.Fluid_ID(i));  mix = string(mp.Mixture(i));
            if ismissing(fid) || fid == "" || any(seen == fid) || any(ismember(split(mix, "-"), cfg.excludeComponents)), continue; end
            seen(end+1) = fid; %#ok<SAGROW>
            zU = str2double(split(erase(string(mp.mole_frac(i)), ["[", "]", " "]), ","));
            zN = str2double(split(erase(string(mp.Nominal_mole_frac(i)), ["[", "]", " "]), ","));
            if numel(zU) ~= numel(zN), continue; end
            hv = split(mix, "-") ~= "C1";
            rel = 100 * max(abs(zU(hv) - zN(hv)) ./ zN(hv));
            aRows(end+1, :) = {fid, mix, string(mp.mole_frac(i)), string(mp.Nominal_mole_frac(i)), rel, ...
                               rel <= 100 * min(cfg.lab.CompRel), rel <= 100 * max(cfg.lab.CompRel)}; %#ok<SAGROW>
        end
        CompCheck = cell2table(aRows, 'VariableNames', {'Fluid_ID', 'Mixture', 'z_used', 'z_nominal', ...
            'HeavyRelChange_pct', 'WithinLow', 'WithinHigh'});
        writetable(CompCheck, xls, 'Sheet', 'CompositionCheck');
        disp(CompCheck);
    end
end

%% 3. BIP windows -------------------------------------------------------------------------
edgeRows = {};  winRows = {};
WIN = struct();                                       % WIN.(pairKey)(q) = [lo, wls, hi, workbook]
for ip = 1:numel(cfg.pairs)
    pr = cfg.pairs{ip};  key = char(strjoin(pr, "_"));  rng = cfg.bipRange.(key);
    sel = find(arrayfun(@(r) all(ismember(pr, split(r.Mixture, "-"))), rows));
    flWB = entities.FluidProperties.loadFromWorkbook(cfg.file, rows(sel(1)).Mixture);
    kWB = flWB.BIPMatrix(flWB.componentIndex(pr(1)), flWB.componentIndex(pr(2)));
    Pfun = @(n, k) localBulk(cfg.file, rock, rows(n), pr, k);
    WIN.(key) = nan(nS, 4);
    for q = 1:nS
        [klo, khi, k1lo, k1hi] = deal(nan(numel(sel), 1));
        for m = 1:numel(sel)
            n = sel(m);
            [lo, hi] = solvers.UncertaintyBudget.band(rows(n).Pm, U(n).UX(q), U(n).UT, U(n).lag, cfg.lab.Repeat_psi * psi);
            klo(m)  = localEdge(@(k) Pfun(n, k) - lo, rng, "lower");
            khi(m)  = localEdge(@(k) Pfun(n, k) - hi, rng, "upper");
            k1lo(m) = localEdge(@(k) Pfun(n, k) - 0.99 * rows(n).Pm, rng, "lower");
            k1hi(m) = localEdge(@(k) Pfun(n, k) - 1.01 * rows(n).Pm, rng, "upper");
            edgeRows(end+1, :) = {strjoin(pr, "-"), 100 * cfg.lab.CompRel(q), rows(n).Label, ...
                (lo - rows(n).Pm) / psi, (hi - rows(n).Pm) / psi, klo(m), khi(m), k1lo(m), k1hi(m)}; %#ok<SAGROW>
        end
        obj = @(k) localWLS(k, sel, rows, U, q, Pfun, cfg.lab.Repeat_psi * psi);
        [kW, chi2] = fminbnd(obj, rng(1), rng(2), optimset('TolX', 1e-6));
        wlo = max(klo);  whi = min(khi);
        WIN.(key)(q, :) = [wlo, kW, whi, kWB];
        status = "ok";
        if any(isnan([klo; khi])), status = "bulk solve failed for a row (see Edges)"; end
        if any(klo == rng(1)) || any(khi == rng(2)), status = status + " | a row does not constrain within the search range"; end
        if wlo > whi, status = "empty window"; end
        winRows(end+1, :) = {strjoin(pr, "-"), 100 * cfg.lab.CompRel(q), numel(sel), wlo, whi, kW, chi2, ...
            kWB, kWB >= wlo && kWB <= whi, max(k1lo), min(k1hi), status}; %#ok<SAGROW>
        fprintf('%-8s %2.0f %%: window [%.4f, %.4f], WLS %.4f (chi2 %.2f), workbook %.4f | +-1 %%: [%.4f, %.4f] %s\n', ...
            strjoin(pr, "-"), 100 * cfg.lab.CompRel(q), wlo, whi, kW, chi2, kWB, max(k1lo), min(k1hi), status);
    end
end
Edges = cell2table(edgeRows, 'VariableNames', {'Pair', 'CompRel_pct', 'Row', 'BandLo_minus_Pm_psi', ...
    'BandHi_minus_Pm_psi', 'BIP_at_BandLo', 'BIP_at_BandHi', 'BIP_at_minus1pct', 'BIP_at_plus1pct'});
Edges = convertvars(Edges, {'Pair', 'Row'}, 'string');
Windows = cell2table(winRows, 'VariableNames', {'Pair', 'CompRel_pct', 'nRows', 'Window_lo', 'Window_hi', ...
    'WLS_BIP', 'WLS_chi2', 'Workbook_BIP', 'Workbook_inside', 'OnePct_lo', 'OnePct_hi', 'Status'});
Windows = convertvars(Windows, {'Pair', 'Status'}, 'string');
writetable(Edges, xls, 'Sheet', 'Edges');
writetable(Windows, xls, 'Sheet', 'Windows');
disp(Windows);

%% 4. Group TIPs at the window levels ----------------------------------------------------------
runner = solvers.DewPointRunner(cfg.file, cfg.poreRadius);
runner.setFWI(fwiRef.k, fwiRef.p, fwiRef.lambda);
groups = thermo.RockOffsetTIP.Groups;
grpPair = ["C1_nC5", "C1_nC8", "C1_nC10", "C1_nC8"];
cache = containers.Map();
tipRows = {};
TIPS = struct();                                       % TIPS.(sKey)(g, level)
for q = find(ismember(cfg.lab.CompRel, cfg.cornerScenarios))
    sKey = sprintf('s%g', 100 * cfg.lab.CompRel(q));
    TIPS.(sKey) = nan(numel(groups), 3);
    for g = 1:numel(groups)
        pk = char(grpPair(g));  w = WIN.(pk)(q, :);
        centre = w(4);  if cfg.cornerCentre == "wls", centre = w(2); end
        lev = [w(1), centre, w(3)];
        if w(1) > w(3) || any(isnan(lev)), lev = [centre, centre, centre]; end   % empty window
        gp = split(groups(g), "|");  pair = split(gp(1), "-");
        idx = dp.findCases(gp(1), gp(2));
        idx = idx(dp.HasConfined(idx) & abs(dp.T(idx) - cfg.T) < 1e-6 & dp.Campaign(idx) ~= "late");
        for L = 1:3
            ck = sprintf('%s|%.8f', groups(g), lev(L));
            if ~isKey(cache, ck)
                runner.clearBIP();  runner.setBIP(pair(1), pair(2), lev(L));
                cache(ck) = localFitGroupTIP(runner, dp, idx, cfg);
            end
            f = cache(ck);
            TIPS.(sKey)(g, L) = f.TIP;
            tipRows(end+1, :) = {100 * cfg.lab.CompRel(q), groups(g), L, lev(L), f.TIP, f.rmse_psi, ...
                strjoin(compose("%+.1f", f.resid_psi.'), " / ")}; %#ok<SAGROW>
            fprintf('  %2.0f %% | %-11s level %d BIP %.5f -> TIP %.4f (RMS %.1f psi)\n', ...
                100 * cfg.lab.CompRel(q), groups(g), L, lev(L), f.TIP, f.rmse_psi);
        end
    end
end
runner.clearBIP();
GroupTIPs = cell2table(tipRows, 'VariableNames', {'CompRel_pct', 'Group', 'Level', 'BIP', 'TIP', ...
    'RMS_resid_psi', 'Residuals_psi'});
writetable(GroupTIPs, xls, 'Sheet', 'GroupTIPs');

%% 5. Corner correlations -------------------------------------------------------------------------
cornRows = {};  rangeRows = {};  CORN = struct();
lvlName = 'LCH';
for q = find(ismember(cfg.lab.CompRel, cfg.cornerScenarios))
    sKey = sprintf('s%g', 100 * cfg.lab.CompRel(q));  Tq = TIPS.(sKey);
    for d = cfg.descriptors
        A = nan(27, 4);  n = 0;  CORN.(sKey).(char(d)) = {};
        for i5 = 1:3, for i8 = 1:3, for i10 = 1:3 %#ok<ALIGN>
            n = n + 1;
            gt = [Tq(1, i5), Tq(2, i8), Tq(3, i10), Tq(4, i8)];
            name = string([lvlName(i5), lvlName(i8), lvlName(i10)]);
            tip = thermo.RockOffsetTIP.custom(name, gt, d, cfg.file, [fwiRef.k, fwiRef.p, fwiRef.lambda]);
            CORN.(sKey).(char(d)){end+1} = struct('name', name, 'lev', [i5, i8, i10], 'tip', tip);
            A(n, :) = [tip.alphaEF2, tip.alphaB1, tip.B, tip.rmseLnEF2];
            cornRows(end+1, :) = {100 * cfg.lab.CompRel(q), d, name, gt(1), gt(2), gt(3), gt(4), ...
                tip.alphaEF2, tip.alphaB1, tip.B, tip.rmseLnEF2}; %#ok<SAGROW>
        end, end, end
        cc = corrcoef(A(:, 1), A(:, 3));
        c = A(14, :);                                  % "CCC" corner
        rangeRows(end+1, :) = {100 * cfg.lab.CompRel(q), d, c(1), min(A(:, 1)), max(A(:, 1)), c(2), ...
            min(A(:, 2)), max(A(:, 2)), c(3), min(A(:, 3)), max(A(:, 3)), cc(1, 2)}; %#ok<SAGROW>
    end
end
Corners = cell2table(cornRows, 'VariableNames', {'CompRel_pct', 'Descriptor', 'Corner', 'TIP_C1nC5_EF2', ...
    'TIP_C1nC8_EF2', 'TIP_C1nC10_EF2', 'TIP_C1nC8_B1', 'alpha_EF2', 'alpha_B1', 'B', 'RMSE_ln_EF2'});
Ranges = cell2table(rangeRows, 'VariableNames', {'CompRel_pct', 'Descriptor', 'alphaEF2_centre', ...
    'alphaEF2_min', 'alphaEF2_max', 'alphaB1_centre', 'alphaB1_min', 'alphaB1_max', 'B_centre', ...
    'B_min', 'B_max', 'corr_alphaEF2_B'});
writetable(Corners, xls, 'Sheet', 'Corners');
writetable(Ranges, xls, 'Sheet', 'CorrelationRanges');
disp(Ranges);

%% 6. Prediction envelope (extreme corners + centre) ------------------------------------------------
if cfg.predict
    psi = entities.DewPointDataset.psi2Pa;
    conf = find(dp.HasConfined & abs(dp.T - cfg.T) < 1e-6 & ~arrayfun(@(i) any(ismember(split(dp.Mixture(i), "-"), cfg.excludeComponents)), ...
                (1:dp.NumCases).'));
    keep = @(lev) all(lev ~= 2) || all(lev == 2);
    predRows = {};
    for q = find(ismember(cfg.lab.CompRel, cfg.cornerScenarios))
        sKey = sprintf('s%g', 100 * cfg.lab.CompRel(q));
        for d = cfg.descriptors
            CC = CORN.(sKey).(char(d));
            for j = 1:numel(CC)
                if ~keep(CC{j}.lev), continue; end
                runner.clearBIP();
                for ip = 1:3
                    pr = cfg.pairs{ip};  w = WIN.(char(strjoin(pr, "_")))(q, :);
                    centre = w(4);  if cfg.cornerCentre == "wls", centre = w(2); end
                    lev = [w(1), centre, w(3)];
                    if w(1) > w(3) || any(isnan(lev)), lev = [centre, centre, centre]; end
                    runner.setBIP(pr(1), pr(2), lev(CC{j}.lev(ip)));
                end
                for i = conf.'
                    c = dp.getCase(i);
                    r = localRadius(c, cfg);
                    fl = entities.FluidProperties.loadFromWorkbook(cfg.file, c.Mixture);
                    K = CC{j}.tip.matrix(fl, c.Rock, c.T);
                    res = runner.solve(c, "FWIOnly", K, 'r', r);
                    Pc = NaN;  if res.ok, Pc = res.P; end
                    predRows(end+1, :) = {100 * cfg.lab.CompRel(q), d, CC{j}.name, c.Mixture, ...
                        "[" + strjoin(compose("%.4g", c.z.'), ",") + "]", c.Rock, r * 1e9, ...
                        c.PdewBulk / psi, res.PdewBulkModel / psi, 100 * (res.PdewBulkModel / c.PdewBulk - 1), ...
                        c.PdewConfMid / psi, c.PdewConfHalfWidth / psi, Pc / psi, 100 * (Pc / c.PdewConfMid - 1), ...
                        res.shift_psi, (c.PdewConfMid - c.PdewBulk) / psi, string(res.reason)}; %#ok<SAGROW>
                end
            end
        end
    end
    runner.clearBIP();
    Pred = cell2table(predRows, 'VariableNames', {'CompRel_pct', 'Descriptor', 'Corner', 'Mixture', 'z', ...
        'Rock', 'rp_nm', 'Bulk_meas_psi', 'Bulk_model_psi', 'Bulk_err_pct', 'Conf_mid_psi', ...
        'Conf_halfwidth_psi', 'Conf_model_psi', 'Conf_err_pct', 'Shift_model_psi', 'Shift_meas_psi', 'Status'});
    writetable(Pred, xls, 'Sheet', 'Predictions');
    % Envelope per case: min / centre (CCC) / max confined error over the corners run
    key = Pred.CompRel_pct + "|" + Pred.Descriptor + "|" + Pred.Mixture + "|" + Pred.z + "|" + Pred.Rock;
    uk = unique(key, 'stable');  envRows = cell(numel(uk), 8);
    for m = 1:numel(uk)
        s_ = key == uk(m);  e_ = Pred.Conf_err_pct(s_);  cn = Pred.Corner(s_);
        ctr = e_(cn == "CCC");  if isempty(ctr), ctr = NaN; end
        i_ = find(s_, 1);
        envRows(m, :) = {Pred.CompRel_pct(i_), Pred.Descriptor(i_), Pred.Mixture(i_), Pred.z(i_), ...
                         Pred.Rock(i_), min(e_), ctr(1), max(e_)};
    end
    env = cell2table(envRows, 'VariableNames', {'CompRel_pct', 'Descriptor', 'Mixture', 'z', 'Rock', ...
        'ConfErr_min_pct', 'ConfErr_centre_pct', 'ConfErr_max_pct'});
    writetable(env, xls, 'Sheet', 'Envelope');
    disp(env);
end

%% 7. Figure: per-row BIP intervals vs the +-1 % rule ---------------------------------------------------
figure('Name', 'BIP windows from the measurement uncertainty', 'Position', [60 60 1400 420]);
tiledlayout(1, numel(cfg.pairs), 'TileSpacing', 'compact');
for ip = 1:numel(cfg.pairs)
    pr = strjoin(cfg.pairs{ip}, "-");  nexttile; hold on; box on;
    rng = cfg.bipRange.(char(strjoin(cfg.pairs{ip}, "_")));
    clip = @(v) min(max(v, rng(1)), rng(2));          % +-Inf edges (no constraint / none) to the range
    E = Edges(Edges.Pair == pr, :);  labs = unique(E.Row, 'stable');
    for m = 1:numel(labs)
        for q = 1:nS
            e = E(E.Row == labs(m) & E.CompRel_pct == 100 * cfg.lab.CompRel(q), :);
            plot(clip([e.BIP_at_BandLo, e.BIP_at_BandHi]), m + [0, 0] + 0.15 * (q - 1), '-', 'LineWidth', 4 - q, ...
                 'DisplayName', sprintf('band, %g %% comp.', 100 * cfg.lab.CompRel(q)), ...
                 'HandleVisibility', localOnOff(m == 1));
        end
        e = E(E.Row == labs(m) & E.CompRel_pct == 100 * cfg.lab.CompRel(1), :);
        plot(clip([e.BIP_at_minus1pct, e.BIP_at_plus1pct]), m - 0.15 + [0, 0], 'k:', 'LineWidth', 2, ...
             'HandleVisibility', localOnOff(m == 1), 'DisplayName', '\pm1 % rule');
    end
    w = Windows(Windows.Pair == pr, :);
    xline(w.Workbook_BIP(1), 'r--', 'workbook', 'HandleVisibility', 'off');
    set(gca, 'YTick', 1:numel(labs), 'YTickLabel', labs, 'TickLabelInterpreter', 'none', 'FontSize', 7);
    xlabel("BIP " + pr);  title(pr);
    if ip == 1, legend('Location', 'best'); end
end
exportgraphics(gcf, fullfile(cfg.outDir, "bip_windows.png"), 'Resolution', 200);
fprintf('\nWritten: %s\n', cfg.outDir);

%% Local functions ------------------------------------------------------------------------------------
function s = localOnOff(tf)
    if tf, s = 'on'; else, s = 'off'; end
end

function P = localBulk(file, rock, row, pair, k)
    % Model bulk Pdew [Pa] of a row with BIP(pair) = k (upper dew branch).
    persistent FL
    if isempty(FL), FL = containers.Map(); end
    key = char(string(file) + "|" + row.Mixture);   % workbook in the key: persists across runs
    if ~isKey(FL, key), FL(key) = entities.FluidProperties.loadFromWorkbook(file, row.Mixture); end
    fl = FL(key);
    K = fl.BIPMatrix;  a = fl.componentIndex(pair(1));  b = fl.componentIndex(pair(2));
    K(a, b) = k;  K(b, a) = k;
    fl = entities.FluidProperties(fl.ComponentNames, fl.Tc, fl.Pc, fl.omega, fl.MW, fl.Vc, ...
        fl.Parachor, fl.ZRackett, fl.VolumeShift, fl.LJ_Size, fl.LJ_Energy, K);
    eos = thermo.ConfinedEOS(fl, rock);
    res = solvers.UncertaintyBudget.upperDew(solvers.SaturationSolver(eos), solvers.PTFlash(eos), ...
                                            row.T, row.Pm, row.z);
    P = NaN;  if res.converged, P = res.P; end
end

function k = localEdge(f, rng, side)
    % BIP at which the model bulk Pdew crosses a band edge (model Pdew rises with
    % the BIP). f = Pmodel - edge. side "lower": BIPs above k satisfy P >= edge;
    % "upper": BIPs below k satisfy P <= edge. No crossing inside rng:
    %   lower: f > 0 everywhere -> rng(1) (no constraint); f < 0 everywhere -> +Inf (none)
    %   upper: f < 0 everywhere -> rng(2) (no constraint); f > 0 everywhere -> -Inf (none)
    % NaN only when the bulk solve fails over the whole scan.
    g = linspace(rng(1), rng(2), 13);  v = arrayfun(f, g);
    for i = 1:numel(g) - 1
        if isfinite(v(i)) && isfinite(v(i + 1)) && v(i) * v(i + 1) < 0
            k = fzero(f, [g(i), g(i + 1)], optimset('TolX', 1e-7));  return;
        end
    end
    v = v(isfinite(v));
    if isempty(v), k = NaN; return; end
    if side == "lower"
        if all(v > 0), k = rng(1); else, k = Inf; end
    else
        if all(v < 0), k = rng(2); else, k = -Inf; end
    end
end

function s = localWLS(k, sel, rows, U, q, Pfun, Urep)
    s = 0;
    for m = 1:numel(sel)
        n = sel(m);
        [t, sg] = solvers.UncertaintyBudget.weight(rows(n).Pm, U(n).UX(q), U(n).UT, U(n).lag, Urep);
        P = Pfun(n, k);
        if ~isfinite(P), s = s + 1e6; else, s = s + ((P - t) / sg)^2; end
    end
end

function r = localRadius(c, cfg)
    r = cfg.poreRadius(char(c.Rock));
    L = cfg.lateCase;
    if c.Mixture == L.Mixture && c.Rock == L.Rock && numel(c.z) == numel(L.z) && max(abs(c.z - L.z)) < 1e-9
        r = cfg.lateRadius_m;
    end
end

function f = localFitGroupTIP(runner, dp, idx, cfg)
    % One group TIP (binary), FWI-Only, k/p/lambda as set on the runner; objective:
    % RMS of the shift error (cfg.tipTarget = "shift") or of the absolute confined error.
    psi = entities.DewPointDataset.psi2Pa;
    resid = @(lt) localGroupResid(runner, dp, idx, -exp(lt), cfg, psi);
    obj = @(lt) sqrt(mean(resid(lt).^2));
    g = cfg.tipSearch.grid;  v = arrayfun(obj, g);
    [~, i] = min(v);
    lo = g(max(i - 1, 1));  hi = g(min(i + 1, numel(g)));
    lt = fminbnd(obj, lo, hi, optimset('TolX', cfg.tipSearch.tolX));
    f = struct('TIP', -exp(lt), 'resid_psi', resid(lt), 'rmse_psi', obj(lt), 'atEdge', i == 1 || i == numel(g));
end

function e = localGroupResid(runner, dp, idx, t, cfg, psi)
    e = nan(numel(idx), 1);
    for n = 1:numel(idx)
        c = dp.getCase(idx(n));
        res = runner.solve(c, "FWIOnly", [0, t; t, 0], 'r', localRadius(c, cfg));
        if ~res.ok, e(n) = 1e3; continue; end            % failure penalty (as in Fit_parameters)
        if cfg.tipTarget == "absolute"
            e(n) = (res.P - c.PdewConfMid) / psi;
        else
            e(n) = res.shift_psi - (c.PdewConfMid - c.PdewBulk) / psi;
        end
    end
end
