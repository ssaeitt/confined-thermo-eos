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
%  The TIP model is NOT hard-coded here: plug the frozen TIP correlation into
%  cfg.kijcFor (section 0). Nothing in this script fits or changes parameters.

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
% Per-case pore-radius override, key "mixture|rock|campaign" (e.g. late B1 core)
cfg.rOverride   = containers.Map();
% cfg.rOverride('C1-nC8|B1|late') = 7.82e-9;

% Frozen FWI parameters (k, p, lambda) of the model being verified
cfg.fwi = struct('k', NaN, 'p', NaN, 'lambda', NaN);

% TIP matrix for a case: @(c, eos, r_cap) -> [nc x nc], symmetric, zero diagonal.
% Plug in the frozen rock-offset correlation used for the HC results.
cfg.kijcFor = @(c, eos, r_cap) error('verify:NoTIP', ...
    'Set cfg.kijcFor to the frozen TIP correlation before running.');

% Bubble-point demonstration feed (liquid-rich; lean feeds have no bubble point)
cfg.bubble = struct('Mixture', "C1-nC10", 'z', [0.5; 0.5], 'P_guess_psi', 2500);

% CCE / CVD case (row index into Mix_Props; [] = first HC case with confined data)
cfg.pvtCaseIndex = [];

% Tolerances
cfg.tol.dP_psi    = 0.01;      % V1/V2 agreement
cfg.tol.lnphi     = 1e-10;     % V4
cfg.flashOffset   = [1e-4, 1e-3];   % V3: P = Pdew*(1 + a), P = Pdew*(1 - b)

if any(isnan([cfg.fwi.k, cfg.fwi.p, cfg.fwi.lambda]))
    error('verify:NoFWI', 'Set cfg.fwi (k, p, lambda) to the frozen model values.');
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
psi = entities.DewPointDataset.psi2Pa;

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
        rows(end+1, :) = {string(label), NaN, NaN, NaN, NaN, NaN, NaN, NaN, "reference failed: " + ref.reason}; %#ok<AGROW>
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
           && eLnPhi < cfg.tol.lnphi;
    note = "";
    if ~s2.converged || abs(dP2) >= cfg.tol.dP_psi
        note = "cold start: " + s2.reason + " (" + s2.route + ")";
    end
    rows(end+1, :) = {string(label), ref.P / psi, ref.shift_psi, dP1, dP2, fa.phases, fb.beta, ...
                      eLnPhi, localVerdict(pass) + " " + note}; %#ok<AGROW>
end
V = cell2table(rows, 'VariableNames', {'Case', 'Pdew_ref_psi', 'shift_psi', 'dP_warm_psi', ...
    'dP_cold_psi', 'phases_above', 'beta_below', 'max_dlnphi_PR', 'Verdict'});
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
    sets = arrayfun(@(e) strjoin(e.ComponentNames, "-"), eqRegions);
    for s = reshape(unique(sets), 1, [])
        grp = eqRegions(sets == s);
        [~, keep] = unique(arrayfun(@(e) e.Label, grp), 'stable');
        grp = grp(keep);
        thermo.EquivalentPR.writeE300(grp, fullfile(expDir, s + "_E300_PROPS.inc"));
        thermo.EquivalentPR.writeGEM(grp,  fullfile(expDir, s + "_GEM_EOS.inc"));
        for e = grp
            e.writeTable(fullfile(expDir, e.Label + ".csv"));
        end
    end
    fprintf('\nSimulator EOS written to %s\n', expDir);
end

%% Local functions ---------------------------------------------------------------------------------
function r = localRadius(cfg, c)
    key = char(c.Mixture + "|" + c.Rock + "|" + c.Campaign);
    if isKey(cfg.rOverride, key)
        r = cfg.rOverride(key);
    else
        r = cfg.poreRadius(char(c.Rock));
    end
end

function s = localVerdict(pass)
    if pass, s = "PASS"; else, s = "FAIL"; end
end
