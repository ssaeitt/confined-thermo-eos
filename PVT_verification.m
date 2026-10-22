%% PVT_VERIFICATION  Verify the general-purpose PVT solvers against the frozen
%  dew-point machinery and export the simulator EOS (FWI-only).
%
%  Frozen model (no placeholders; values as in TIP_correlation.m / Baclcalc_rp_B1_late.m):
%    preset jL (primary)    : joint-Linearized k, p, lambda; TIP form LL-inv_wbar
%    preset w1 (sensitivity): pareto-exact-w1 k, p, lambda; TIP form LL-Tstar
%    TIP_ij = -exp(alpha_rock + B * x_ij), rock-offset, size-independent.
%    alpha_EF2, alpha_B1, B are re-derived here by the SAME least-squares fit as
%    TIP_correlation.m (fitForm) from the frozen per-group optimised TIPs, using the
%    component properties of the workbook. No TIP workbook is needed.
%
%  Checks, per HC case with confined data (CO2 mixtures excluded), at every pore
%  radius assigned to the case (late B1 core: nominal AND back-calculated):
%   V1  SaturationSolver dew point == DewPointRunner (FlashEngine) dew point
%   V2  SaturationSolver dew point from a cold start == V1
%   V3  PTFlash: single phase just above Pdew, two phases just below
%   V4  EquivalentPR reproduces ln(phi) and Z of ConfinedEOS to round-off
%  Then:
%   B   bulk and confined bubble points of liquid-rich feeds (+ flash check).
%       The SIGN of the confined bubble-point shift is a model prediction (the
%       model is calibrated on dew points only); it is reported, never imposed.
%   C   CCE / CVD, bulk vs confined, for one case
%   E   simulator EOS export (ECLIPSE 300 / CMG GEM templates + CSV)
%
%  Needs (local, not in the repository): config/MixtureData - Complete Mixtures.xlsx
%  Nothing in this script fits model parameters to data.

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
cfg.mixtureFile = fullfile(scriptDir, "config", "MixtureData - Complete Mixtures.xlsx");
cfg.outDir      = fullfile(scriptDir, "results", "pvt_verification");
cfg.T           = 293.15;
cfg.kB          = 1.380649e-23;                                % [J/K]
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});   % nominal radii [m]

% --- Frozen model ---------------------------------------------------------------
cfg.preset = "jL";                     % "jL" (primary) | "w1" (sensitivity)
presets.jL = struct('combination', "joint-linearized", ...
    'k', 112.42634826677, 'p', 1.3558179290678, 'lambda', 2.0916235090532, ...
    'descriptor', "inv_wbar", 'form', "LL-inv_wbar", ...
    'tipGroups', [-4.9229696302796, -1.96920366521686, -1.44628069292089, -0.754862049795788]);
presets.w1 = struct('combination', "pareto-exact-w1", ...
    'k', 166.410386, 'p', 1.548842, 'lambda', 4.479082, ...    % as in TIP_correlation.m cfg.combos
    'descriptor', "Tstar", 'form', "LL-Tstar", ...
    'tipGroups', [-7.012375, -2.883946, -2.124570, -1.168401]);
% Groups of the per-group optimised TIPs, in the order of tipGroups (FWI-Only)
cfg.tipGroups = ["C1-nC5|EF2", "C1-nC8|EF2", "C1-nC10|EF2", "C1-nC8|B1"];
ps = presets.(cfg.preset);

% Late B1 core, C1-nC8 [0.9, 0.1]: nominal AND back-calculated radii (band edges kept
% for reference). Values from co2_pipeline_v3.m (rpManual, jL).
cfg.lateCase = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.lateRadii_nm = struct('nominal', 11.25, 'backcalc', 7.82711826617637, ...
                          'backcalc_lo', 7.62188369515433, 'backcalc_hi', 8.04796163088063);
cfg.lateRadiiToRun = ["nominal", "backcalc"];     % add "backcalc_lo", "backcalc_hi" if wanted

% Bubble-point demonstration feeds (liquid-rich; lean feeds have no bubble point)
cfg.bubble = struct('Mixture', {"C1-nC10", "C1-nC8"}, 'z', {[0.5; 0.5], [0.5; 0.5]}, ...
                    'P_guess_psi', {2500, 2500});

% CCE / CVD case (row index into Mix_Props; [] = first HC case with confined data)
cfg.pvtCaseIndex = [];

% Tolerances
cfg.tol.dP_psi  = 0.01;            % V1/V2 agreement
cfg.tol.round   = 1e-10;           % V4 (ln phi and Z)
cfg.flashOffset = [1e-4, 1e-3];    % V3: P = Pdew*(1 + a), P = Pdew*(1 - b)

if ~isfile(cfg.mixtureFile)
    error('verify:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.mixtureFile);
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

%% 1. Data, engines and the rock-offset TIP correlation ------------------------------
dp = entities.DewPointDataset.loadFromWorkbook(cfg.mixtureFile, 'DefaultT', cfg.T);
isHC = ~contains(dp.Mixture, "CO2");
idx  = find(dp.HasConfined & isHC).';
runner = solvers.DewPointRunner(cfg.mixtureFile, cfg.poreRadius);
for i = idx, runner.engine(dp.Mixture(i), dp.Rock(i)); end
runner.setFWI(ps.k, ps.p, ps.lambda);

tip = localFitTIP(cfg, ps, runner);
fprintf('Model %s (%s): k = %.6f, p = %.6f, lambda = %.6f\n', cfg.preset, ps.combination, ps.k, ps.p, ps.lambda);
fprintf('TIP %s: alpha_EF2 = %+.6f, alpha_B1 = %+.6f, B = %+.6f (refit RMSE ln, EF2 = %.4f)\n\n', ...
        ps.form, tip.alphaEF2, tip.alphaB1, tip.B, tip.rmseLnEF2);
kijcFor = @(c, eos) localTIPMatrix(eos.Fluid, c.Rock, c.T, ps.descriptor, tip, cfg.kB);

%% 2. V1-V4 per case and radius ------------------------------------------------------
rows = {};
eqRegions = thermo.EquivalentPR.empty(1, 0);
for i = idx
    c   = dp.getCase(i);
    eng = runner.engine(c.Mixture, c.Rock);
    K   = kijcFor(c, eng.eos);
    [radii, rLabels] = localRadii(cfg, c);

    for q = 1:numel(radii)
        r = radii(q);
        label = sprintf('%s | %s | z1 = %.4f | r_p %s', c.Mixture, c.Rock, c.z(1), rLabels(q));
        ref = runner.solve(c, "FWIOnly", K, 'r', r, 'Key', string(label));   % sets eng.eos.kijc
        eos = eng.eos;
        if ~ref.ok
            rows(end+1, :) = {string(label), r * 1e9, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                              "reference failed: " + ref.reason}; %#ok<AGROW>
            continue;
        end

        stab = solvers.PhaseStability(eos);
        sat  = solvers.SaturationSolver(eos, 'Stability', stab);
        pt   = solvers.PTFlash(eos, 'Stability', stab);

        s1 = sat.solve("dew", c.T, ref.P, c.z, r, 'K_seed', exp(ref.lnK));   % V1 warm
        s2 = sat.solve("dew", c.T, 1.05 * c.PdewBulk, c.z, r);               % V2 cold
        fa = pt.solve(ref.P * (1 + cfg.flashOffset(1)), c.T, c.z, r);        % V3
        fb = pt.solve(ref.P * (1 - cfg.flashOffset(2)), c.T, c.z, r, 'K_seed', exp(ref.lnK));

        eq = thermo.EquivalentPR.fromEOS(eos, c.T, r, sprintf("%s-%s-r%.2fnm", c.Rock, c.Mixture, r * 1e9));
        Pgrid = linspace(0.3, 1.5, 7).' * ref.P;
        zgrid = [c.z, fb.x, fb.y, ones(numel(c.z), 1) / numel(c.z)];
        [eLnPhi, eZ] = eq.verify(eos, Pgrid, zgrid);                         % V4
        eqRegions(end+1) = eq; %#ok<AGROW>

        dP1 = (s1.P - ref.P) / psi;
        dP2 = (s2.P - ref.P) / psi;
        pass = s1.converged && abs(dP1) < cfg.tol.dP_psi && fa.phases == 1 && fb.phases == 2 ...
               && eLnPhi < cfg.tol.round && eZ < cfg.tol.round;
        note = "";
        if ~s2.converged || abs(dP2) >= cfg.tol.dP_psi
            note = "cold start: " + s2.reason + " (" + s2.route + ")";
        end
        rows(end+1, :) = {string(label), r * 1e9, ref.P / psi, ref.shift_psi, dP1, dP2, fa.phases, ...
                          fb.beta, eLnPhi, eZ, localVerdict(pass) + " " + note}; %#ok<AGROW>
    end
end
V = cell2table(rows, 'VariableNames', {'Case', 'rp_nm', 'Pdew_ref_psi', 'shift_psi', 'dP_warm_psi', ...
    'dP_cold_psi', 'phases_above', 'beta_below', 'max_dlnphi_PR', 'max_dZ_PR', 'Verdict'});
disp(V);
writetable(V, fullfile(cfg.outDir, "verification_" + cfg.preset + ".csv"));

%% 3. B: bubble points of liquid-rich feeds -------------------------------------------------
fprintf('\nBubble points at %.2f K (sign of the confined shift is a model prediction)\n', cfg.T);
bubRows = {};
for b = cfg.bubble
    for rock = ["EF2", "B1"]
        cb  = struct('Mixture', b.Mixture, 'Rock', rock, 'z', b.z, 'T', cfg.T);
        eng = runner.engine(cb.Mixture, cb.Rock);
        r   = cfg.poreRadius(char(rock));
        eng.eos.kijc = kijcFor(cb, eng.eos);
        eng.eos.ConfinementScale = 1;
        sat = solvers.SaturationSolver(eng.eos);
        pt  = solvers.PTFlash(eng.eos);
        bb  = sat.solve("bubble", cfg.T, b.P_guess_psi * psi, cb.z, Inf);
        bc  = sat.solve("bubble", cfg.T, bb.P, cb.z, r, 'K_seed', bb.K);
        if ~(bb.converged && bc.converged)
            fprintf('  %-8s %-4s bulk: %s | confined: %s\n', b.Mixture, rock, bb.reason, bc.reason);
            bubRows(end+1, :) = {b.Mixture, rock, r * 1e9, NaN, NaN, NaN, NaN, NaN, ...
                                 "bulk: " + bb.reason + " | confined: " + bc.reason}; %#ok<AGROW>
            continue;
        end
        fa = pt.solve(bc.P * 1.001, cfg.T, cb.z, r);
        fb = pt.solve(bc.P * 0.999, cfg.T, cb.z, r, 'K_seed', bc.K);
        fprintf('  %-8s %-4s Pb bulk %8.2f | confined %8.2f psi (shift %+8.2f) | flash +0.1%%: %d ph, -0.1%%: %d ph\n', ...
                b.Mixture, rock, bb.P / psi, bc.P / psi, (bc.P - bb.P) / psi, fa.phases, fb.phases);
        bubRows(end+1, :) = {b.Mixture, rock, r * 1e9, bb.P / psi, bc.P / psi, (bc.P - bb.P) / psi, ...
                             fa.phases, fb.phases, "ok"}; %#ok<AGROW>
    end
end
if ~isempty(bubRows)
    Bt = cell2table(bubRows, 'VariableNames', {'Mixture', 'Rock', 'rp_nm', 'Pb_bulk_psi', ...
        'Pb_conf_psi', 'shift_psi', 'phases_above', 'phases_below', 'status'});
    writetable(Bt, fullfile(cfg.outDir, "bubble_points_" + cfg.preset + ".csv"));
end

%% 4. C: CCE / CVD, bulk vs confined ----------------------------------------------------------
pvtIdx = cfg.pvtCaseIndex;
if isempty(pvtIdx), pvtIdx = idx(1); end
c   = dp.getCase(pvtIdx);
[radii, ~] = localRadii(cfg, c);
r   = radii(1);
eng = runner.engine(c.Mixture, c.Rock);
ref = runner.solve(c, "FWIOnly", kijcFor(c, eng.eos), 'r', r);
if ref.ok
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
else
    warning('verify:PVTCase', 'CCE/CVD skipped: reference dew point failed (%s).', ref.reason);
end

%% 5. E: simulator EOS export ------------------------------------------------------------------
% One deck per component set (all EOS regions in a deck must share components).
expDir = fullfile(cfg.outDir, "simulator_eos_" + cfg.preset);
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

%% Local functions ------------------------------------------------------------------------------
function [radii, labels] = localRadii(cfg, c)
    % Late B1 core: every radius in cfg.lateRadiiToRun; other cases: nominal rock radius.
    L = cfg.lateCase;
    isLate = c.Mixture == L.Mixture && c.Rock == L.Rock && numel(c.z) == numel(L.z) ...
             && max(abs(c.z(:) - L.z(:))) < 1e-9;
    if isLate
        labels = cfg.lateRadiiToRun;
        radii  = zeros(1, numel(labels));
        for k = 1:numel(labels)
            radii(k) = cfg.lateRadii_nm.(char(labels(k))) * 1e-9;
        end
    else
        labels = "nominal";
        radii  = cfg.poreRadius(char(c.Rock));
    end
end

function tip = localFitTIP(cfg, ps, runner)
    % Same fit as TIP_correlation.m (fitForm): ln(-TIP) = alpha_EF2*[EF2] + alpha_B1*[B1] + B*x
    G = numel(cfg.tipGroups);
    x = zeros(G, 1);  isB1 = false(G, 1);
    for g = 1:G
        pr  = split(cfg.tipGroups(g), "|");
        eng = runner.engine(pr(1), pr(2));
        X   = localDescriptor(eng.eos.Fluid, ps.descriptor, cfg.T, cfg.kB);
        x(g) = X(1, 2);
        isB1(g) = pr(2) == "B1";
    end
    y = log(-ps.tipGroups(:));
    A = [double(~isB1), double(isB1), x];
    if rank(A) < 3, error('verify:TIPRank', 'TIP design matrix is rank-deficient.'); end
    par = A \ y;
    res = A * par - y;
    tip = struct('alphaEF2', par(1), 'alphaB1', par(2), 'B', par(3), ...
                 'rmseLnEF2', sqrt(mean(res(~isB1).^2)), 'xCal', x);
end

function K = localTIPMatrix(fl, rock, T, descriptor, tip, kB)
    switch string(rock)
        case "EF2", alpha = tip.alphaEF2;
        case "B1",  alpha = tip.alphaB1;
        otherwise,  error('verify:Rock', 'No rock offset for %s.', rock);
    end
    nc = fl.NC;
    K = -exp(alpha + tip.B * localDescriptor(fl, descriptor, T, kB));
    K(1:nc+1:end) = 0;
    K = (K + K.') / 2;
end

function X = localDescriptor(fl, name, T, kB)
    % Pair descriptors, identical to TIP_correlation.m pairDescriptor (LJ_Energy in J).
    w  = fl.omega(:);
    e  = fl.LJ_Energy(:);
    MW = fl.MW(:);
    Ts = (kB * T) ./ sqrt(e * e.');
    switch name
        case "wbar",      X = 0.5 * (w + w.');
        case "inv_wbar",  X = 2 ./ (w + w.');
        case "Tstar",     X = Ts;
        case "inv_Tstar", X = 1 ./ Ts;
        case "ln_Tstar",  X = log(Ts);
        case "inv_Mgeo",  X = 1 ./ sqrt(MW * MW.');
        otherwise, error('verify:Descriptor', 'Unknown descriptor "%s".', name);
    end
end

function s = localVerdict(pass)
    if pass, s = "PASS"; else, s = "FAIL"; end
end
