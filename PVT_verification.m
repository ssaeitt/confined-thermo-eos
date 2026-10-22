%% PVT_VERIFICATION  Verify the general-purpose PVT solvers against the frozen
%  dew-point machinery and export the simulator EOS (FWI-only).
%
%  Frozen model: thermo.RockOffsetTIP (single source of k, p, lambda, the per-group
%  TIPs and the rock-offset correlation fit; presets jL = primary, w1 = sensitivity).
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
%   C   CCE / CVD, bulk vs confined, for the cases in cfg.pvtCases; dropout also on
%       the laboratory basis (vol% of V_sat) with MLDO vs the measured bulk MLDO
%   E   simulator EOS export (ECLIPSE 300 / CMG GEM templates + CSV), including a
%       bulk region (r_p = Inf) per component set for fractures / dead volume
%
%  Needs (local, not in the repository): config/MixtureData - Complete Mixtures.xlsx
%  Nothing in this script fits model parameters to data.

clear; clc; close all;

%% Preflight ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% 0. Configuration -------------------------------------------------------------
% --- Model variant -----------------------------------------------------------------
%   "frozen" : current calibration (workbook "MixtureData - Complete Mixtures.xlsx", preset jL)
%   "recal"  : C1-nC5 recalibration (workbook "MixtureData - Recalibrated.xlsx", preset jL2)
cfg.variant = "frozen";
variants = struct( ...
    'frozen', struct('file', "MixtureData - Complete Mixtures.xlsx", 'preset', "jL"), ...
    'recal',  struct('file', "MixtureData - Recalibrated.xlsx",     'preset', "jL2"));
vsel = variants.(char(cfg.variant));
cfg.mixtureFile = fullfile(scriptDir, "config", vsel.file);
cfg.outDir      = fullfile(scriptDir, "results", "pvt_verification_" + cfg.variant);
cfg.T           = 293.15;
cfg.kB          = 1.380649e-23;                                % [J/K]
cfg.poreRadius  = containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9});   % nominal radii [m]

% --- Frozen model ---------------------------------------------------------------
cfg.preset = vsel.preset;              % "jL" | "jL2" (from cfg.variant) | "w1" (sensitivity)

% Late B1 core, C1-nC8 [0.9, 0.1]: nominal AND back-calculated radii (band edges kept
% for reference). Values from co2_pipeline_v3.m (rpManual, jL).
cfg.lateCase = struct('Mixture', "C1-nC8", 'Rock', "B1", 'z', [0.9; 0.1]);
cfg.lateRadii_nm = struct('nominal', 11.25, 'backcalc', 7.82711826617637, ...
                          'backcalc_lo', 7.62188369515433, 'backcalc_hi', 8.04796163088063);
cfg.lateRadiiToRun = ["nominal", "backcalc"];     % add "backcalc_lo", "backcalc_hi" if wanted

% Bubble-point demonstration feeds (liquid-rich; lean feeds have no bubble point)
cfg.bubble = struct('Mixture', {"C1-nC10", "C1-nC8"}, 'z', {[0.5; 0.5], [0.5; 0.5]}, ...
                    'P_guess_psi', {2500, 2500});

% CCE / CVD cases: rows of Mix_Props by Fluid_ID and rock (late B1 core: both radii)
cfg.pvtCases = struct('Fluid', {"F2", "F7", "F8", "F5"}, 'Rock', {"EF2", "B1", "B1", "EF2"}, ...
                      'Campaign', {"early", "early", "late", "early"});
cfg.PminFrac = 0.10;               % CCE/CVD down to 10 % of Psat

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
mpRaw = readtable(cfg.mixtureFile, 'Sheet', 'Mix_Props', 'VariableNamingRule', 'preserve');
hasMLDO = all(ismember(["Fluid_ID", "Campaign", "MLDO_pct", "P_MLDO_psi"], string(mpRaw.Properties.VariableNames)));

tip = thermo.RockOffsetTIP.frozen(cfg.preset, cfg.mixtureFile);
fprintf('%s\n\n', tip.summary());
runner.setFWI(tip.k, tip.p, tip.lambda);
kijcFor = @(c, eos) tip.matrix(eos.Fluid, c.Rock, c.T);

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
% Dropout on the laboratory basis (vol% of V_sat); MLDO compared with the measured BULK
% MLDO (confined MLDO has not been measured). Psat values are the model dew points.
if ~hasMLDO
    warning('verify:NoMLDO', 'Mix_Props has no Fluid_ID/MLDO columns: CCE/CVD skipped.');
    cfg.pvtCases = cfg.pvtCases([]);
end
mldoRows = {};
for pc = cfg.pvtCases
    sel = find(string(mpRaw.Fluid_ID) == pc.Fluid & string(mpRaw.Rock) == pc.Rock & ...
               string(mpRaw.Campaign) == pc.Campaign, 1);
    if isempty(sel)
        warning('verify:PVTCase', 'No Mix_Props row for %s on %s (%s).', pc.Fluid, pc.Rock, pc.Campaign);
        continue;
    end
    c = dp.getCase(sel);                         % Mix_Props and dp share row order
    [radii, rLabels] = localRadii(cfg, c);
    eng = runner.engine(c.Mixture, c.Rock);
    pvt = solvers.PVTExperiments(solvers.PTFlash(eng.eos));
    stem = sprintf('%s_%s', pc.Fluid, pc.Rock);
    if pc.Campaign == "late", stem = stem + "_late"; end

    curves = {};
    for q = 1:numel(radii)
        ref = runner.solve(c, "FWIOnly", kijcFor(c, eng.eos), 'r', radii(q));
        if ~ref.ok
            warning('verify:PVTCase', '%s r_p %s: reference dew point failed (%s).', stem, rLabels(q), ref.reason);
            continue;
        end
        if q == 1                                % bulk once per case
            Pb = ref.PdewBulkModel;
            Pl = solvers.PVTExperiments.pressureGrid(Pb, cfg.PminFrac * Pb, 40);
            cceB = pvt.cce(c.T, c.z, Inf, Pb, [1.05 * Pb; Pl]);
            [cvdB, ~] = pvt.cvd(c.T, c.z, Inf, Pb, Pl);
            writetable(cceB, fullfile(cfg.outDir, "cce_" + stem + "_bulk.csv"));
            writetable(cvdB, fullfile(cfg.outDir, "cvd_" + stem + "_bulk.csv"));
            [m, pm] = solvers.PVTExperiments.mldo(cceB);
            mldoRows(end+1, :) = {pc.Fluid, pc.Rock, pc.Campaign, "bulk", Inf, Pb / psi, m, pm, ...
                                  mpRaw.MLDO_pct(sel), mpRaw.P_MLDO_psi(sel)}; %#ok<AGROW>
            curves(end+1, :) = {"bulk", cceB, cvdB}; %#ok<AGROW>
        end
        Pl = solvers.PVTExperiments.pressureGrid(ref.P, cfg.PminFrac * ref.P, 40);
        cceC = pvt.cce(c.T, c.z, radii(q), ref.P, [1.05 * ref.P; Pl]);
        [cvdC, ~] = pvt.cvd(c.T, c.z, radii(q), ref.P, Pl);
        tagR = sprintf('r%.2fnm', radii(q) * 1e9);
        writetable(cceC, fullfile(cfg.outDir, "cce_" + stem + "_" + tagR + ".csv"));
        writetable(cvdC, fullfile(cfg.outDir, "cvd_" + stem + "_" + tagR + ".csv"));
        [m, pm] = solvers.PVTExperiments.mldo(cceC);
        mldoRows(end+1, :) = {pc.Fluid, pc.Rock, pc.Campaign, "confined " + rLabels(q), radii(q) * 1e9, ...
                              ref.P / psi, m, pm, mpRaw.MLDO_pct(sel), mpRaw.P_MLDO_psi(sel)}; %#ok<AGROW>
        curves(end+1, :) = {sprintf('confined, r_p = %.2f nm', radii(q) * 1e9), cceC, cvdC}; %#ok<AGROW>
    end
    if isempty(curves), continue; end

    figure('Name', "Liquid dropout " + stem, 'Position', [60 60 1150 430]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; hold on; box on;
    for k = 1:size(curves, 1)
        plot(curves{k, 2}.P_psi, curves{k, 2}.LDO_pct_Vsat, '.-', 'DisplayName', curves{k, 1});
    end
    if isfinite(mpRaw.MLDO_pct(sel)) && isfinite(mpRaw.P_MLDO_psi(sel))
        plot(mpRaw.P_MLDO_psi(sel), mpRaw.MLDO_pct(sel), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
             'DisplayName', 'measured bulk MLDO');
    end
    xlabel('P [psi]'); ylabel('CCE liquid dropout [% of V_{sat}]'); legend('Location', 'best');
    title(sprintf('%s (%s, z_1 = %.4f), %s', pc.Fluid, c.Mixture, c.z(1), c.Rock));
    nexttile; hold on; box on;
    for k = 1:size(curves, 1)
        plot(curves{k, 3}.P_psi, 100 * curves{k, 3}.SL_cell, '.-', 'DisplayName', curves{k, 1});
    end
    xlabel('P [psi]'); ylabel('CVD liquid dropout [% of V_{cell}]'); legend('Location', 'best');
    exportgraphics(gcf, fullfile(cfg.outDir, "dropout_" + stem + ".png"), 'Resolution', 200);
end
if ~isempty(mldoRows)
    M = cell2table(mldoRows, 'VariableNames', {'Fluid', 'Rock', 'Campaign', 'State', 'rp_nm', ...
        'Psat_model_psi', 'MLDO_model_pct', 'P_MLDO_model_psi', 'MLDO_exp_bulk_pct', 'P_MLDO_exp_bulk_psi'});
    M = convertvars(M, {'Fluid', 'Rock', 'Campaign', 'State'}, 'string');
    disp(M);
    writetable(M, fullfile(cfg.outDir, "mldo_" + cfg.preset + ".csv"));
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
        % Bulk region (r_p = Inf): standard PR with the workbook BIPs, listed first
        eosB = runner.engine(strjoin(grp(1).ComponentNames, "-"), cfg.lateCase.Rock).eos;
        grp  = [thermo.EquivalentPR.fromEOS(eosB, cfg.T, Inf, "bulk-" + setName), grp]; %#ok<AGROW>
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

function s = localVerdict(pass)
    if pass, s = "PASS"; else, s = "FAIL"; end
end
