%% PVT_VERIFICATION  Verify the new general-purpose PVT solvers against the
%  frozen dew-point machinery, and export the simulator EOS (FWI-only).
%
%  Checks, per HC case with confined data (CO2 mixtures are excluded):
%   V1  SaturationSolver dew point == DewPointRunner (FlashEngine) dew point
%   V2  SaturationSolver dew point from a cold start (homotopy) == V1
%   V3  PTFlash: single phase just above Pdew, two phases just below
%   V4  EquivalentPR reproduces ln(phi) and Z of ConfinedEOS to round-off
%  Then:
%   B   bulk and confined bubble points of a liquid-rich feed + flash check
%   C   CCE / CVD, bulk vs confined, for one case
%   E   simulator EOS export (ECLIPSE 300 / CMG GEM templates + CSV)
%
%  Model (frozen, FWI-Only), same definitions as Baclcalc_rp_B1_late.m:
%    k, p, lambda from cfg.preset (jL = joint-Linearized, primary)
%    TIP_ij = -exp(alpha_rock + B * x_ij), rock-offset, size-independent,
%    coefficients read from config/tip_correlation_rockoffset.xlsx
%    (sheet "Correlation"), primary form LL-inv_wbar, x_ij = 2/(omega_i + omega_j).
%  Every binary TIP is gated against the "Pair_TIPs" sheet before use.
%  Nothing in this script fits or changes parameters.

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
cfg.mixtureFile = fullfile(scriptDir, "config", "MixtureData - Complete Mixtures.xlsx");
cfg.outDir      = fullfile(scriptDir, "results", "pvt_verification");
cfg.T           = 293.15;
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});
% Per-case pore-radius override, key "mixture|rock|z1,z2,..." (exact composition)
cfg.rOverride   = containers.Map();

% Frozen model preset (values as in Baclcalc_rp_B1_late.m)
cfg.preset = "jL";
presets.jL = struct('combination', "joint-linearized", 'k', 112.42634826677, 'p', 1.3558179290678, ...
    'lambda', 2.0916235090532, 'form', "LL-inv_wbar");
presets.w1 = struct('combination', "pareto-exact-w1", 'k', 166.410385934507, 'p', 1.5488421885245, ...
    'lambda', 4.47908205631447, 'form', "LL-Tstar");
ps = presets.(cfg.preset);
cfg.fwi = struct('k', ps.k, 'p', ps.p, 'lambda', ps.lambda);
cfg.tip = struct('file', fullfile(scriptDir, "config", "tip_correlation_rockoffset.xlsx"), ...
    'combination', ps.combination, 'mode', "FWIOnly", 'form', ps.form, 'gateTol', 1e-6);
cfg.kB = 1.380649e-23;                                   % [J/K]

% Late B1 core: C1-nC8 [0.9, 0.1] uses the back-calculated radius (7.82 nm).
cfg.rOverride('C1-nC8|B1|0.9,0.1') = 7.82e-9;

% Bubble-point demonstration feed (liquid-rich; lean feeds have no bubble point)
cfg.bubble = struct('Mixture', "C1-nC10", 'z', [0.5; 0.5], 'P_guess_psi', 2500);

% CCE / CVD case (row index into Mix_Props; [] = first HC case with confined data)
cfg.pvtCaseIndex = [];

% Tolerances
cfg.tol.dP_psi    = 0.01;      % V1/V2 agreement
cfg.tol.lnphi     = 1e-10;     % V4
cfg.flashOffset   = [1e-4, 1e-3];   % V3: P = Pdew*(1 + a), P = Pdew*(1 - b)

if ~isfile(cfg.tip.file)
    error('verify:NoTIPFile', 'TIP correlation workbook not found: %s', cfg.tip.file);
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

tipTab  = readtable(cfg.tip.file, 'Sheet', 'Correlation', 'VariableNamingRule', 'preserve');
pairTab = readtable(cfg.tip.file, 'Sheet', 'Pair_TIPs',   'VariableNamingRule', 'preserve');
cfg.kijcFor = @(c, eos, r_cap) localTIP(c, eos, cfg, tipTab, pairTab);
fprintf('Model: %s | k = %.6f, p = %.6f, lambda = %.6f | TIP %s (%s)\n\n', cfg.preset, ...
        cfg.fwi.k, cfg.fwi.p, cfg.fwi.lambda, cfg.tip.form, cfg.tip.combination);

%% 1. Data and engines -------------------------------------------------------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
isHC = ~contains(dp.Mixture, "CO2");
idx  = find(dp.HasConfined & isHC).';
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius);
for i = idx, runner.engine(dp.Mixture(i), dp.Rock(i)); end
runner.setFWI(cfg.fwi.k, cfg.fwi.p, cfg.fwi.lambda);

rOf = @(c) localRadius(cfg, c);

%% 2. V1-V4 per case ------------------------------------------------------------------
rows = {};
eqRegions = thermo.EquivalentPR.empty(1, 0);
for i = idx
    c   = dp.getCase(i);
    r   = rOf(c);
    eng = runner.engine(c.Mixture, c.Rock);
    K   = cfg.kijcFor(c, eng.eos, r);
    ref = runner.solve(c, "FWIOnly", K, 'r', r);         % sets eng.eos.kijc
    eos = eng.eos;

    stab = solvers.PhaseStability(eos);
    sat  = solvers.SaturationSolver(eos, 'Stability', stab);
    pt   = solvers.PTFlash(eos, 'Stability', stab);

    label = sprintf('%s | %s | z1 = %.4f', c.Mixture, c.Rock, c.z(1));
    if ~ref.ok
        rows(end+1, :) = {string(label), r * 1e9, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "reference failed: " + ref.reason}; %#ok<AGROW>
        continue;
    end

    % V1: warm-started from the reference K
    s1 = sat.solve("dew", c.T, ref.P, c.z, r, 'K_seed', exp(ref.lnK));
    % V2: cold start (TPD/Wilson seed, homotopy fallback)
    s2 = sat.solve("dew", c.T, 1.05 * c.PdewBulk, c.z, r);

    % V3: flash just above and just below the confined dew point
    fa = pt.solve(ref.P * (1 + cfg.flashOffset(1)), c.T, c.z, r);
    fb = pt.solve(ref.P * (1 - cfg.flashOffset(2)), c.T, c.z, r, 'K_seed', exp(ref.lnK));

    % V4: exact PR equivalent
    eq = thermo.EquivalentPR.fromEOS(eos, c.T, r, sprintf("%s-%s-r%.2fnm", c.Rock, c.Mixture, r * 1e9));
    Pgrid = linspace(0.3, 1.5, 7).' * ref.P;
    zgrid = [c.z, fb.x, fb.y, ones(numel(c.z), 1) / numel(c.z)];
    [eLnPhi, eZ] = eq.verify(eos, Pgrid, zgrid);
    eqRegions(end+1) = eq; %#ok<AGROW>

    dP1 = (s1.P - ref.P) / psi;
    dP2 = (s2.P - ref.P) / psi;
    pass = s1.converged && abs(dP1) < cfg.tol.dP_psi && fa.phases == 1 && fb.phases == 2 ...
           && eLnPhi < cfg.tol.lnphi && eZ < cfg.tol.lnphi;
    note = "";
    if ~s2.converged || abs(dP2) >= cfg.tol.dP_psi
        note = "cold start: " + s2.reason + " (" + s2.route + ")";
    end
    rows(end+1, :) = {string(label), r * 1e9, ref.P / psi, ref.shift_psi, dP1, dP2, fa.phases, fb.beta, ...
                      eLnPhi, eZ, localVerdict(pass) + " " + note}; %#ok<AGROW>
end
V = cell2table(rows, 'VariableNames', {'Case', 'rp_nm', 'Pdew_ref_psi', 'shift_psi', 'dP_warm_psi', ...
    'dP_cold_psi', 'phases_above', 'beta_below', 'max_dlnphi_PR', 'max_dZ_PR', 'Verdict'});
disp(V);
writetable(V, fullfile(cfg.outDir, "verification.csv"));

%% 3. B: bubble points of a liquid-rich feed ----------------------------------------------
fprintf('\nBubble points, %s, z = [%s]\n', cfg.bubble.Mixture, num2str(cfg.bubble.z.', '%.3f '));
for rock = ["EF2", "B1"]
    cb  = struct('Mixture', cfg.bubble.Mixture, 'Rock', rock, 'z', cfg.bubble.z, 'T', cfg.T, ...
                 'Campaign', "", 'PdewBulk', NaN);
    eng = runner.engine(cb.Mixture, cb.Rock);
    r   = cfg.poreRadius(char(rock));
    eng.eos.kijc = cfg.kijcFor(cb, eng.eos, r);
    eng.eos.ConfinementScale = 1;
    sat = solvers.SaturationSolver(eng.eos);
    pt  = solvers.PTFlash(eng.eos);
    bb  = sat.solve("bubble", cfg.T, cfg.bubble.P_guess_psi * psi, cb.z, Inf);
    bc  = sat.solve("bubble", cfg.T, bb.P, cb.z, r, 'K_seed', bb.K);
    if ~(bb.converged && bc.converged)
        fprintf('  %-4s bulk: %s | confined: %s\n', rock, bb.reason, bc.reason);
        continue;
    end
    fa = pt.solve(bc.P * 1.001, cfg.T, cb.z, r);
    fb = pt.solve(bc.P * 0.999, cfg.T, cb.z, r, 'K_seed', bc.K);
    fprintf('  %-4s Pb bulk %8.2f psi | confined %8.2f psi (shift %+7.2f) | flash +0.1%%: %d ph, -0.1%%: %d ph, beta = %.2e\n', ...
            rock, bb.P / psi, bc.P / psi, (bc.P - bb.P) / psi, fa.phases, fb.phases, fb.beta);
end

%% 4. C: CCE / CVD, bulk vs confined ---------------------------------------------------------
pvtIdx = cfg.pvtCaseIndex;
if isempty(pvtIdx), pvtIdx = idx(1); end
c   = dp.getCase(pvtIdx);
r   = rOf(c);
eng = runner.engine(c.Mixture, c.Rock);
K   = cfg.kijcFor(c, eng.eos, r);
ref = runner.solve(c, "FWIOnly", K, 'r', r);
pt  = solvers.PTFlash(eng.eos);
pvt = solvers.PVTExperiments(pt);
Pb  = ref.PdewBulkModel;
Plist_b = linspace(0.98, 0.15, 25).' * Pb;
Plist_c = linspace(0.98, 0.15, 25).' * ref.P;
[cvdB, ~] = pvt.cvd(c.T, c.z, Inf, Pb, Plist_b);
[cvdC, ~] = pvt.cvd(c.T, c.z, r,   ref.P, Plist_c);
cceB = pvt.cce(c.T, c.z, Inf, Pb, [1.2 * Pb; Plist_b]);
cceC = pvt.cce(c.T, c.z, r,   ref.P, [1.2 * ref.P; Plist_c]);
writetable(cvdB, fullfile(cfg.outDir, "cvd_bulk.csv"));
writetable(cvdC, fullfile(cfg.outDir, "cvd_confined.csv"));
writetable(cceB, fullfile(cfg.outDir, "cce_bulk.csv"));
writetable(cceC, fullfile(cfg.outDir, "cce_confined.csv"));

figure('Name', 'Liquid dropout');
tiledlayout(1, 2);
nexttile; hold on; box on;
plot(cceB.P_psi, 100 * cceB.SL, 'o-', 'DisplayName', 'bulk');
plot(cceC.P_psi, 100 * cceC.SL, 's-', 'DisplayName', sprintf('confined, r_p = %.2f nm', r * 1e9));
xlabel('P [psi]'); ylabel('Liquid dropout, CCE [% of V_{total}]'); legend('Location', 'best');
title(sprintf('%s, %s, z_1 = %.3f', c.Mixture, c.Rock, c.z(1)));
nexttile; hold on; box on;
plot(cvdB.P_psi, 100 * cvdB.SL_cell, 'o-', 'DisplayName', 'bulk');
plot(cvdC.P_psi, 100 * cvdC.SL_cell, 's-', 'DisplayName', 'confined');
xlabel('P [psi]'); ylabel('Liquid dropout, CVD [% of V_{cell}]'); legend('Location', 'best');
exportgraphics(gcf, fullfile(cfg.outDir, "dropout_bulk_vs_confined.png"), 'Resolution', 200);

%% 5. E: simulator EOS export ----------------------------------------------------------------------
% One deck per component set (all EOS regions in a deck must share components).
expDir = fullfile(cfg.outDir, "simulator_eos");
if ~isfolder(expDir), mkdir(expDir); end
if ~isempty(eqRegions)
    sets   = strings(1, numel(eqRegions));
    labels = strings(1, numel(eqRegions));
    for k = 1:numel(eqRegions)
        sets(k)   = strjoin(eqRegions(k).ComponentNames, "-");
        labels(k) = eqRegions(k).Label;
    end
    for setName = unique(sets, 'stable')
        inSet = find(sets == setName);
        [~, keep] = unique(labels(inSet), 'stable');      % one region per rock and r_p
        grp = eqRegions(inSet(keep));
        thermo.EquivalentPR.writeE300(grp, fullfile(expDir, setName + "_E300_PROPS.inc"));
        thermo.EquivalentPR.writeGEM(grp,  fullfile(expDir, setName + "_GEM_EOS.inc"));
        for k = 1:numel(grp)
            grp(k).writeTable(fullfile(expDir, grp(k).Label + ".csv"));
        end
    end
    fprintf('\nSimulator EOS written to %s\n', expDir);
end

%% Local functions ---------------------------------------------------------------------------------
function r = localRadius(cfg, c)
    % Priority: explicit override > r_p column of Mix_Props > nominal rock radius.
    key = char(c.Mixture + "|" + c.Rock + "|" + strjoin(string(c.z(:).'), ","));
    if isKey(cfg.rOverride, key)
        r = cfg.rOverride(key);
    elseif isfield(c, 'PoreRadius_m') && isfinite(c.PoreRadius_m)
        r = c.PoreRadius_m;
    else
        r = cfg.poreRadius(char(c.Rock));
    end
end

function s = localVerdict(pass)
    if pass, s = "PASS"; else, s = "FAIL"; end
end

function K = localTIP(c, eos, cfg, tipTab, pairTab)
    % Rock-offset TIP matrix, ln(-TIP_ij) = alpha_rock + B * x_ij, for every pair
    % of the case's components; binary pairs gated against the Pair_TIPs sheet.
    sel = string(tipTab.Combination) == cfg.tip.combination & ...
          string(tipTab.Mode) == cfg.tip.mode & string(tipTab.Form) == cfg.tip.form;
    if nnz(sel) ~= 1
        error('verify:CorrelationRow', '%d Correlation rows match %s / %s / %s.', nnz(sel), ...
              cfg.tip.combination, cfg.tip.mode, cfg.tip.form);
    end
    vn = "alpha_" + string(c.Rock);
    if ~ismember(vn, string(tipTab.Properties.VariableNames))
        error('verify:NoRockOffset', 'Column %s not found in the Correlation sheet.', vn);
    end
    alpha = tipTab.(vn)(sel);
    B     = tipTab.B(sel);
    descr = string(tipTab.Descriptor(sel));

    fl = eos.Fluid;
    w  = fl.omega(:);
    e  = fl.LJ_Energy(:);
    MW = fl.MW(:);
    Ts = (cfg.kB * c.T) ./ sqrt(e * e.');
    switch descr
        case "wbar",      X = 0.5 * (w + w.');
        case "inv_wbar",  X = 2 ./ (w + w.');
        case "Tstar",     X = Ts;
        case "inv_Tstar", X = 1 ./ Ts;
        case "ln_Tstar",  X = log(Ts);
        case "inv_Mgeo",  X = 1 ./ sqrt(MW * MW.');
        otherwise
            error('verify:UnknownDescriptor', 'Descriptor "%s" is not implemented.', descr);
    end
    nc = fl.NC;
    K = -exp(alpha + B * X);
    K(1:nc+1:end) = 0;
    K = (K + K.') / 2;

    % Gate: binary pairs must reproduce the Pair_TIPs sheet for this rock
    if nc == 2
        nm  = string(fl.ComponentNames);
        sel = string(pairTab.Combination) == cfg.tip.combination & string(pairTab.Mode) == cfg.tip.mode & ...
              string(pairTab.Form) == cfg.tip.form & string(pairTab.Context) == "binary" & ...
              string(pairTab.Rock) == string(c.Rock) & string(pairTab.Pair) == nm(1) + "-" + nm(2);
        if nnz(sel) == 1
            ref = pairTab.TIP(sel);
            if abs(K(1, 2) - ref) > cfg.tip.gateTol * max(1, abs(ref))
                error('verify:TIPGate', 'TIP(%s, %s-%s) = %.8f vs Pair_TIPs %.8f.', ...
                      c.Rock, nm(1), nm(2), K(1, 2), ref);
            end
        end
    end
end
