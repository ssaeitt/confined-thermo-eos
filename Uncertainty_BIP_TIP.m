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
% 1b. GUM / NIST TN 1297 form (cfg.gum): every bound becomes a standard
%    uncertainty (rectangular: bound/sqrt(3); one-sided lag and dew pick:
%    correction +w/2, u = w/sqrt(12)), combined in quadrature, expanded with
%    k = 2 (~0.95 level of confidence). Reported as y +- U_c with JCED digits,
%    with E_n = (model - y)/U_c (|E_n| <= 1: model consistent with the datum).
%    Sheets GUM_Pdew, GUM_MLDO, GUM_MLDO_dT, GUM_Confined, GUM_Duplicates, GUM_Footnote.
%    Duplicates (same mixture and z, separate batches) give a Type A check:
%    |dP| / (k sqrt(2) u_c) <= 1 if the Type B budget covers the repeatability.
%    cfg.gum.windows = true builds the BIP windows and WLS weights from y +- U_c.
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
cfg.lab = struct('CompRel', [0.01, 0.05, 0.10], 'TempRel', 0.05, 'dT_C', 1, 'CellVolume_cc', 20, ...
                 'DetectVolume_cc', 0.1, 'Repeat_psi', 0, 'Gauge_psi', 0, 'GaugeRel', 0, ...
                 'PickWidth_psi', 0, 'SameGauge', true);
% Laboratory answers (27 Oct 2026): temperature +-1 degC both at room temperature
% (5 % of 20 degC) and in the oven (40, 60 degC) -> dT_C = 1 (bound, rectangular);
% composition "below 1 %" (gas mixture) against the earlier 5-10 % -> all three kept
% as scenarios; transducer reading 4015-4020 for 4000 psi -> not entered until checked
% (it may be the 14.7 psi absolute/gauge offset rather than a span error).
% GaugeRel      : transducer span bound as a fraction of the reading (e.g. 0.004)
% SameGauge     : bulk (cell) and confined (core) read on the same transducer, so its
%                 offset cancels in the shift and only the span acts (on the shift)
% dT_C          : temperature figure in K; overrides TempRel when finite (e.g. 0.1 from
%                 the thermostat/RTD calibration, with cfg.gum.TempDist = "standard")
% Gauge_psi     : pressure transducer accuracy (bound unless GaugeDist = "standard")
% PickWidth_psi : one-sided dew-pick width = depletion rate x observation interval
%                 (e.g. 1 psi/min read every 5 min -> 5 psi); 0 = not quantified
% GUM settings (sheets GUM_*). Distributions: "rectangular" | "triangular" | "normal95" | "standard"
cfg.gum = struct('enable', true, 'Coverage', 2, 'CompDist', "rectangular", 'TempDist', "rectangular", ...
                 'ReadDist', "rectangular", 'GaugeDist', "rectangular", ...
                 'repeat_psi', 0, ...                % Type A standard uncertainty added to every row
                 'repeatFromDuplicates', false, ...  % true: use the pooled s of the duplicates instead
                 'windows', true);                   % BIP windows / WLS from y +- U_c (false: bound bands)
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
cfg.cornerScenarios = [0.01, 0.10];               % composition scenarios for the corners
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
ub = solvers.UncertaintyBudget('CompRel', cfg.lab.CompRel, 'TempRel', cfg.lab.TempRel, 'dT_C', cfg.lab.dT_C, ...
        'CellVolume_cc', cfg.lab.CellVolume_cc, 'DetectVolume_cc', cfg.lab.DetectVolume_cc, ...
        'CompDist', cfg.gum.CompDist, 'TempDist', cfg.gum.TempDist, 'ReadDist', cfg.gum.ReadDist, ...
        'GaugeDist', cfg.gum.GaugeDist, 'Gauge_psi', cfg.lab.Gauge_psi, 'GaugeRel', cfg.lab.GaugeRel, ...
        'PickWidth_psi', cfg.lab.PickWidth_psi, 'Coverage', cfg.gum.Coverage);
gumOn = cfg.gum.enable;
if ~gumOn, cfg.gum.windows = false; end
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
vn = {'Case', 'Pmeas_psi', 'Pmodel_psi', 'Bulk_err_pct', 'dP_Tm1_psi', 'dP_Tp1_psi', 'UT_psi'};
for q = 1:nS
    tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
    vn = [vn, {['UX' tag '_psi'], ['U' tag '_psi'], ['U' tag '_pct'], ['UX' tag '_oneSided']}]; %#ok<AGROW>
end
vn = [vn, {'lag_psi', 'lag_pct', 'Note'}];
Budget = cell2table(budRows, 'VariableNames', vn);
Budget = convertvars(Budget, {'Case', 'Note'}, 'string');
writetable(Budget, xls, 'Sheet', 'Budget');
disp(Budget);

%% 2a. GUM form: duplicates (Type A check), y +- U_c, E_n -----------------------------------------
G = [];
if gumOn
    % Duplicates: same mixture and z, separate fillings (rows differ in measured Pdew)
    keyR = strings(nR, 1);
    for n = 1:nR, keyR(n) = rows(n).Mixture + "|" + strjoin(compose("%.6g", rows(n).z.'), ","); end
    [uKey, ~, gi] = unique(keyR, 'stable');
    dupRows = {};  ssq = 0;  dof = 0;
    for gk = 1:numel(uKey)
        mem = find(gi == gk);
        if numel(mem) < 2, continue; end
        Pv = [rows(mem).Pm].';  sg = std(Pv);              % Type A, n - 1
        ssq = ssq + (numel(mem) - 1) * sg^2;  dof = dof + numel(mem) - 1;
        g1 = ub.gum(U(mem(1)).P0, rows(mem(1)).Pm, U(mem(1)).UX, U(mem(1)).UT, U(mem(1)).lag, 0);
        dPd = max(Pv) - min(Pv);
        rec = {uKey(gk), numel(mem), strjoin(compose("%.0f", Pv.' / psi), " / "), dPd / psi, sg / psi, numel(mem) - 1};
        for q = 1:nS
            rec = [rec, {g1.uTypeB(q) / psi, dPd / (ub.Coverage * sqrt(2) * g1.uTypeB(q))}]; %#ok<AGROW>
        end
        dupRows(end+1, :) = rec; %#ok<SAGROW>
        fprintf('  duplicate %-30s dP %.0f psi, s %.1f psi (nu %d) | u_c(Type B) %s psi | dP/(k sqrt2 u_c) %s\n', ...
            uKey(gk), dPd / psi, sg / psi, numel(mem) - 1, strjoin(compose("%.1f", g1.uTypeB / psi), " / "), ...
            strjoin(compose("%.2f", dPd ./ (ub.Coverage * sqrt(2) * g1.uTypeB)), " / "));
    end
    sPool = NaN;  if dof > 0, sPool = sqrt(ssq / dof); end
    uRep = cfg.gum.repeat_psi * psi;
    if cfg.gum.repeatFromDuplicates && isfinite(sPool), uRep = sPool; end
    if ~isempty(dupRows)
        vnD = {'Group', 'n', 'Pmeas_psi', 'Range_psi', 'sTypeA_psi', 'nu'};
        for q = 1:nS
            tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
            vnD = [vnD, {['ucTypeB' tag '_psi'], ['Ratio' tag]}]; %#ok<AGROW>
        end
        Dups = cell2table(dupRows, 'VariableNames', vnD);
        Dups = convertvars(Dups, {'Group', 'Pmeas_psi'}, 'string');
        writetable(Dups, xls, 'Sheet', 'GUM_Duplicates');
        fprintf('  pooled Type A s = %.1f psi (nu = %d); u_rep used = %.1f psi\n', sPool / psi, dof, uRep / psi);
    end

    G = repmat(struct('y', NaN, 'uc', nan(1, nS), 'Uc', nan(1, nS), 'En', nan(1, nS)), nR, 1);
    gRows = {};  ucr = ub.uCompRel();
    for n = 1:nR
        g = ub.gum(U(n).P0, rows(n).Pm, U(n).UX, U(n).UT, U(n).lag, uRep);
        G(n) = struct('y', g.y, 'uc', g.uc, 'Uc', g.Uc, 'En', g.En);
        hv = split(rows(n).Mixture, "-") ~= "C1";
        rec = {rows(n).Label, rows(n).Mixture, "[" + strjoin(compose("%.4g", rows(n).z.'), ",") + "]", ...
               rows(n).T, ub.uTemp(rows(n).T), rows(n).Pm / psi, g.y / psi, U(n).P0 / psi, g.dev / psi, ...
               g.uT / psi, g.uG / psi, g.uLag / psi, g.uPick / psi, g.uRep / psi};
        for q = 1:nS
            rec = [rec, {ucr(q) * max(rows(n).z(hv)), g.uX(q) / psi, g.uc(q) / psi, g.Uc(q) / psi, ...
                         100 * g.Uc(q) / g.y, g.En(q), ...
                         solvers.UncertaintyBudget.report(g.y / psi, g.Uc(q) / psi), ...
                         solvers.UncertaintyBudget.report(g.y / psi * 1e-3 * 6.894757, g.Uc(q) / psi * 1e-3 * 6.894757)}]; %#ok<AGROW>
        end
        gRows(end+1, :) = rec; %#ok<SAGROW>
    end
    vnG = {'Case', 'Mixture', 'z', 'T_K', 'uT_K', 'Pmeas_psi', 'y_psi', 'Pmodel_psi', 'Dev_psi', ...
           'uT_psi', 'uGauge_psi', 'uLag_psi', 'uPick_psi', 'uRep_psi'};
    for q = 1:nS
        tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
        vnG = [vnG, {['ux_heavy' tag], ['uX' tag '_psi'], ['uc' tag '_psi'], ['Uc' tag '_psi'], ['Uc' tag '_pct'], ...
                     ['En' tag], ['Report' tag '_psi'], ['Report' tag '_MPa']}]; %#ok<AGROW>
    end
    GumP = cell2table(gRows, 'VariableNames', vnG);
    GumP = convertvars(GumP, {'Case', 'Mixture', 'z'}, 'string');
    writetable(GumP, xls, 'Sheet', 'GUM_Pdew');
    disp(GumP(:, [{'Case', 'y_psi', 'Pmodel_psi'}, vnG(contains(vnG, ["Uc", "En"]) & ~contains(vnG, "pct"))]));

    % Footnote in the JCED/TRC style (values at the first row's temperature)
    Tf = rows(1).T;
    ft = sprintf(['Standard uncertainties: u(T) = %.2g K (%s, from %s); u_r(x_heavy) = %s (%s, from %s %% ' ...
        'relative on the heavy component as charged); u(p_gauge) = %.2g psi + %.2g %% of reading; detection lag (%.2g cm3 in %.3g cm3) and ' ...
        'dew pick (%.3g psi) one-sided: value corrected by +w/2, u = w/sqrt(12); u_rep = %.2g psi. ' ...
        'Expanded uncertainty U_c(p_dew) = k u_c with k = %g (%s level of confidence).'], ...
        ub.uTemp(Tf), cfg.gum.TempDist, localTempText(cfg, Tf), strjoin(compose("%.3g", ub.uCompRel), " / "), ...
        cfg.gum.CompDist, strjoin(compose("%g", 100 * cfg.lab.CompRel), " / "), ...
        cfg.lab.Gauge_psi / solvers.UncertaintyBudget.divisor(cfg.gum.GaugeDist), ...
        100 * cfg.lab.GaugeRel / solvers.UncertaintyBudget.divisor(cfg.gum.GaugeDist), cfg.lab.DetectVolume_cc, ...
        cfg.lab.CellVolume_cc, cfg.lab.PickWidth_psi, uRep / psi, ub.Coverage, localConf(ub.Coverage));
    writetable(table(string(ft), 'VariableNames', {'Footnote'}), xls, 'Sheet', 'GUM_Footnote');
    fprintf('\n%s\n\n', ft);
end

%% 2b. MLDO budget (Bulk_PVT: all fluids and temperatures) ----------------------------------
if cfg.mldo
    bp = readtable(cfg.file, 'Sheet', 'Bulk_PVT', 'VariableNamingRule', 'preserve');
    mRows = {};  gmRows = {};
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
        if gumOn
            gm = ub.gumMLDO(m.M0, bp.MLDO_pct(i), m.UX, m.UT);
            grec = {string(bp.Fluid_ID(i)), mix, bp.T(i), ub.uTemp(bp.T(i)), bp.MLDO_pct(i), m.M0, gm.dev, ...
                    gm.uT, gm.uRead};
            for q = 1:nS
                grec = [grec, {gm.uX(q), gm.uc(q), gm.Uc(q), gm.En(q), ...
                               solvers.UncertaintyBudget.report(gm.y, gm.Uc(q))}]; %#ok<AGROW>
            end
            gmRows(end+1, :) = grec; %#ok<SAGROW>
        end
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
    if gumOn && ~isempty(gmRows)
        vnGM = {'Fluid_ID', 'Mixture', 'T_K', 'uT_K', 'MLDO_meas_pct', 'MLDO_model_pct', 'Dev_pts', 'uT_pts', 'uRead_pts'};
        for q = 1:nS
            tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
            vnGM = [vnGM, {['uX' tag '_pts'], ['uc' tag '_pts'], ['Uc' tag '_pts'], ['En' tag], ['Report' tag]}]; %#ok<AGROW>
        end
        GumM = cell2table(gmRows, 'VariableNames', vnGM);
        GumM = convertvars(GumM, {'Fluid_ID', 'Mixture'}, 'string');
        writetable(GumM, xls, 'Sheet', 'GUM_MLDO');
        disp(GumM(:, [{'Fluid_ID', 'T_K', 'MLDO_meas_pct', 'MLDO_model_pct'}, vnGM(startsWith(vnGM, ["Uc", "En"]))]));

        % Temperature differences of one fluid (same batch): the composition term
        % cancels to first order; u^2 = uT(T)^2 + uT(Tref)^2 + 2 uRead^2.
        dRows = {};
        for fid = unique(GumM.Fluid_ID, 'stable').'
            sF = find(GumM.Fluid_ID == fid);
            if numel(sF) < 2, continue; end
            [~, o] = min(GumM.T_K(sF));  r0 = sF(o);
            for r1 = setdiff(sF, r0).'
                dMe = GumM.MLDO_meas_pct(r1) - GumM.MLDO_meas_pct(r0);
                dMo = GumM.MLDO_model_pct(r1) - GumM.MLDO_model_pct(r0);
                uD  = sqrt(GumM.uT_pts(r1)^2 + GumM.uT_pts(r0)^2 + GumM.uRead_pts(r1)^2 + GumM.uRead_pts(r0)^2);
                dRows(end+1, :) = {fid, GumM.T_K(r0), GumM.T_K(r1), dMe, dMo, dMo - dMe, ub.Coverage * uD, ...
                                   (dMo - dMe) / (ub.Coverage * uD), ...
                                   solvers.UncertaintyBudget.report(dMe, ub.Coverage * uD)}; %#ok<SAGROW>
            end
        end
        if ~isempty(dRows)
            GumMdT = cell2table(dRows, 'VariableNames', {'Fluid_ID', 'Tref_K', 'T_K', 'dMLDO_meas_pts', ...
                'dMLDO_model_pts', 'Dev_pts', 'Uc_pts', 'En', 'Report'});
            writetable(GumMdT, xls, 'Sheet', 'GUM_MLDO_dT');
            disp(GumMdT);
        end
    end

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
            if cfg.gum.windows
                lo = G(n).y - G(n).Uc(q);  hi = G(n).y + G(n).Uc(q);
            else
                [lo, hi] = solvers.UncertaintyBudget.band(rows(n).Pm, U(n).UX(q), U(n).UT, U(n).lag, cfg.lab.Repeat_psi * psi);
            end
            klo(m)  = localEdge(@(k) Pfun(n, k) - lo, rng, "lower");
            khi(m)  = localEdge(@(k) Pfun(n, k) - hi, rng, "upper");
            k1lo(m) = localEdge(@(k) Pfun(n, k) - 0.99 * rows(n).Pm, rng, "lower");
            k1hi(m) = localEdge(@(k) Pfun(n, k) - 1.01 * rows(n).Pm, rng, "upper");
            edgeRows(end+1, :) = {strjoin(pr, "-"), 100 * cfg.lab.CompRel(q), rows(n).Label, ...
                (lo - rows(n).Pm) / psi, (hi - rows(n).Pm) / psi, klo(m), khi(m), k1lo(m), k1hi(m)}; %#ok<SAGROW>
        end
        Gw = [];  if cfg.gum.windows, Gw = G; end
        obj = @(k) localWLS(k, sel, rows, U, q, Pfun, cfg.lab.Repeat_psi * psi, Gw);
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
Edges = cell2table(edgeRows, 'VariableNames', {'Pair', 'CompRel_pct', 'Case', 'BandLo_minus_Pm_psi', ...
    'BandHi_minus_Pm_psi', 'BIP_at_BandLo', 'BIP_at_BandHi', 'BIP_at_minus1pct', 'BIP_at_plus1pct'});
Edges = convertvars(Edges, {'Pair', 'Case'}, 'string');
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
        centre = w(4);  if cfg.cornerCentre == "wls" || w(4) < w(1) || w(4) > w(3), centre = w(2); end   % workbook outside window -> WLS
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

%% 4b. GUM form of the confined data (preset model, workbook BIPs) -------------------------------
% Measured range [lower, upper] read as rectangular: mid, u_range = half-width / sqrt(3).
% Absolute confined Pdew: u^2 = u_range^2 + u_X^2 + u_T^2 + u_gauge^2 + u_rep^2, with the
%   composition and temperature sensitivities of the bulk row (shift << Pdew).
% Shift (confined - bulk): same batch in the cell and the core, so the composition
%   term cancels to first order; temperature, gauge and repeatability enter twice
%   (two apparatus), the cell detection lag and dew pick once (bulk side only).
if gumOn
    confIdx = find(dp.HasConfined & abs(dp.T - cfg.T) < 1e-6 & ...
        ~arrayfun(@(i) any(ismember(split(dp.Mixture(i), "-"), cfg.excludeComponents)), (1:dp.NumCases).'));
    cRows = {};
    for i = confIdx.'
        c = dp.getCase(i);
        n = find(arrayfun(@(r) r.Mixture == c.Mixture && numel(r.z) == numel(c.z) && ...
                 max(abs(r.z - c.z)) < 1e-9 && abs(r.Pm - c.PdewBulk) < 1e-6, rows), 1);
        if isempty(n), continue; end
        g = ub.gum(U(n).P0, rows(n).Pm, U(n).UX, U(n).UT, U(n).lag, uRep);
        uR = c.PdewConfHalfWidth / sqrt(3);
        fl = entities.FluidProperties.loadFromWorkbook(cfg.file, c.Mixture);
        res = runner.solve(c, "FWIOnly", fwiRef.matrix(fl, c.Rock, c.T), 'r', localRadius(c, cfg));
        Pc = NaN;  Sc = NaN;  if res.ok, Pc = res.P;  Sc = res.shift_psi * psi; end
        Smeas = c.PdewConfMid - g.y;
        if cfg.lab.SameGauge          % offset cancels; span acts on the shift only
            uGs2 = (cfg.lab.GaugeRel * abs(Smeas) / solvers.UncertaintyBudget.divisor(cfg.gum.GaugeDist))^2;
        else
            uGs2 = 2 * g.uG^2;
        end
        uS = sqrt(uR^2 + 2 * g.uT^2 + uGs2 + 2 * g.uRep^2 + g.uLag^2 + g.uPick^2);
        rec = {c.Mixture, "[" + strjoin(compose("%.4g", c.z.'), ",") + "]", c.Rock, c.Campaign, ...
               localRadius(c, cfg) * 1e9, c.PdewConfMid / psi, uR / psi, Pc / psi, Smeas / psi, Sc / psi, ...
               ub.Coverage * uS / psi, (Sc - Smeas) / (ub.Coverage * uS), ...
               solvers.UncertaintyBudget.report(Smeas / psi, ub.Coverage * uS / psi)};
        for q = 1:nS
            uA = sqrt(uR^2 + g.uX(q)^2 + g.uT^2 + g.uG^2 + g.uRep^2);
            rec = [rec, {ub.Coverage * uA / psi, (Pc - c.PdewConfMid) / (ub.Coverage * uA), ...
                         solvers.UncertaintyBudget.report(c.PdewConfMid / psi, ub.Coverage * uA / psi)}]; %#ok<AGROW>
        end
        cRows(end+1, :) = rec; %#ok<SAGROW>
    end
    vnC = {'Mixture', 'z', 'Rock', 'Campaign', 'rp_nm', 'Conf_mid_psi', 'uRange_psi', 'Conf_model_psi', ...
           'Shift_meas_psi', 'Shift_model_psi', 'Uc_shift_psi', 'En_shift', 'Report_shift_psi'};
    for q = 1:nS
        tag = sprintf('%g', 100 * cfg.lab.CompRel(q));
        vnC = [vnC, {['Uc' tag '_abs_psi'], ['En' tag '_abs'], ['Report' tag '_abs_psi']}]; %#ok<AGROW>
    end
    if ~isempty(cRows)
        GumC = cell2table(cRows, 'VariableNames', vnC);
        GumC = convertvars(GumC, {'Mixture', 'z', 'Rock', 'Campaign'}, 'string');
        writetable(GumC, xls, 'Sheet', 'GUM_Confined');
        disp(GumC(:, [{'Mixture', 'z', 'Rock', 'Shift_meas_psi', 'Shift_model_psi', 'Uc_shift_psi', 'En_shift'}, ...
                      vnC(startsWith(vnC, "En") & endsWith(vnC, "_abs"))]));
    end
end

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
                    centre = w(4);  if cfg.cornerCentre == "wls" || w(4) < w(1) || w(4) > w(3), centre = w(2); end   % workbook outside window -> WLS
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
bandCol = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];   % fixed colour per composition scenario
for ip = 1:numel(cfg.pairs)
    pr = strjoin(cfg.pairs{ip}, "-");  nexttile; hold on; box on;
    rng = cfg.bipRange.(char(strjoin(cfg.pairs{ip}, "_")));
    clip = @(v) min(max(v, rng(1)), rng(2));          % +-Inf edges (no constraint / none) to the range
    E = Edges(Edges.Pair == pr, :);  labs = unique(E.Case, 'stable');
    for m = 1:numel(labs)
        for q = 1:nS
            e = E(E.Case == labs(m) & E.CompRel_pct == 100 * cfg.lab.CompRel(q), :);
            plot(clip([e.BIP_at_BandLo, e.BIP_at_BandHi]), m + [0, 0] + 0.15 * (q - 1), '-', 'LineWidth', 4 - q, ...
                 'Color', bandCol(min(q, size(bandCol, 1)), :), ...
                 'DisplayName', sprintf('band, %g %% comp.', 100 * cfg.lab.CompRel(q)), ...
                 'HandleVisibility', localOnOff(m == 1));
        end
        e = E(E.Case == labs(m) & E.CompRel_pct == 100 * cfg.lab.CompRel(1), :);
        plot(clip([e.BIP_at_minus1pct, e.BIP_at_plus1pct]), m - 0.15 + [0, 0], 'k:', 'LineWidth', 2, ...
             'HandleVisibility', localOnOff(m == 1), 'DisplayName', '\pm1 % rule');
    end
    w = Windows(Windows.Pair == pr, :);
    xline(w.Workbook_BIP(1), 'r--', 'workbook', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom', ...
          'LabelOrientation', 'horizontal', 'FontSize', 7);
    ylim([0.5, numel(labs) + 0.6]);
    set(gca, 'YTick', 1:numel(labs), 'YTickLabel', labs, 'TickLabelInterpreter', 'none', 'FontSize', 7);
    xlabel("BIP " + pr);  title(pr);
    if ip == 1, legend('Location', 'best'); end
end
exportgraphics(gcf, fullfile(cfg.outDir, "bip_windows.png"), 'Resolution', 200);
fprintf('\nWritten: %s\n', cfg.outDir);

%% Local functions ------------------------------------------------------------------------------------
function s = localTempText(cfg, T)
    if isfinite(cfg.lab.dT_C), s = sprintf('%g K', cfg.lab.dT_C);
    else, s = sprintf('%g %% of the reading in degC = %.2g K at %.2f K', 100 * cfg.lab.TempRel, ...
                      cfg.lab.TempRel * (T - 273.15), T); end
end

function s = localConf(k)
    % Level of confidence for a normal distribution and coverage factor k.
    s = sprintf('%.3g', erf(k / sqrt(2)));
end

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

function s = localWLS(k, sel, rows, U, q, Pfun, Urep, G)
    % G non-empty: GUM target y and standard uncertainty u_c; else the bound form.
    s = 0;
    for m = 1:numel(sel)
        n = sel(m);
        if ~isempty(G)
            t = G(n).y;  sg = G(n).uc(q);
        else
            [t, sg] = solvers.UncertaintyBudget.weight(rows(n).Pm, U(n).UX(q), U(n).UT, U(n).lag, Urep);
        end
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
