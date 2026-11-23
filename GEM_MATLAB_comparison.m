%% GEM_MATLAB_COMPARISON  CCE of C1-nC8 (z = 0.923 / 0.077) at the pressures of
%  the GEM no-flow blocks, for bulk, EF2 54 nm and B1 11.25 nm (jL3, FWI-only).
%
%  Reads  tools/gem_comparison/pressures_all.csv       (one column, P_psia: the
%         pressures of the 225 columns of GEM deck tools/gem_comparison/CCEV_C_jL3.dat)
%  Writes results/gem_comparison/matlab_<tag>.csv      one row per block pressure:
%         P_psia, phases, SL, rhoL_kgm3, rhoV_kgm3, x_nC8, y_nC8, Z_L, Z_V, converged
%         results/gem_comparison/matlab_dewpoints.csv  dew point per set
%  The same pressures are the GEM block pressures, so GEM, MATLAB and the Python
%  reference (tools/gem_comparison/python_<tag>.csv) can be compared
%  with tools/gem_comparison/read_gem_cce.py and compare_gem_matlab.py.
%
%  SL = liquid volume / total volume at P, the same basis as GEM's oil saturation
%  in a block with no water. Densities are mass densities with the volume shift.
%  Nothing here fits parameters: the frozen jL3 model is used as calibrated.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir, "src"));

%% Configuration ---------------------------------------------------------------
cfg.mixtureFile = fullfile(scriptDir, "config", "MixtureData - Recalibrated psia.xlsx");
cfg.preset  = "jL3";
cfg.T       = 293.15;                      % 68 F, as in the GEM decks
cfg.z       = [0.923; 0.077];              % C1, nC8
cfg.inFile  = fullfile(scriptDir, "tools", "gem_comparison", "pressures_all.csv");
cfg.dir     = fullfile(scriptDir, "results", "gem_comparison");
cfg.Pguess_psi = 4100;                     % starting guess for every dew point
sets = struct('tag',  {"bulk", "EF2_54", "B1_11.25"}, ...
              'Rock', {"EF2",  "EF2",    "B1"}, ...
              'r',    {Inf,    54e-9,    11.25e-9});

psi = entities.DewPointDataset.psi2Pa;
if ~isfile(cfg.mixtureFile)
    error('cmp:NoWorkbook', 'Workbook not found (config/ is local only): %s', cfg.mixtureFile);
end
if ~isfolder(cfg.dir), mkdir(cfg.dir); end
Pall = sort(readtable(cfg.inFile).P_psia, 'descend');   % march down so K can be reused

%% Engines and the frozen rock term ---------------------------------------------
runner = solvers.DewPointRunner(cfg.mixtureFile, containers.Map({'EF2', 'B1'}, {54e-9, 11.25e-9}));
tip = thermo.RockOffsetTIP.frozen(cfg.preset, cfg.mixtureFile);
runner.setFWI(tip.k, tip.p, tip.lambda);

dew = table('Size', [0 3], 'VariableTypes', {'string', 'double', 'string'}, ...
            'VariableNames', {'set', 'Pdew_psia', 'route'});

for s = sets
    eng = runner.engine("C1-nC8", s.Rock);
    eos = eng.eos;
    tip.applyTo(eos);                     % frozen k, p, lambda (as set on the runner)
    % Rock term for this rock (has no effect when r = Inf, where the fluid-wall term is zero).
    eos.kijc = tip.matrix(eos.Fluid, s.Rock, cfg.T);

    stab = solvers.PhaseStability(eos);
    sat  = solvers.SaturationSolver(eos, 'Stability', stab);
    pt   = solvers.PTFlash(eos, 'Stability', stab);

    d = sat.solve("dew", cfg.T, cfg.Pguess_psi * psi, cfg.z, s.r);
    if ~d.converged, warning('cmp:Dew', '%s: dew point did not converge (%s)', s.tag, d.reason); end
    dew(end+1, :) = {s.tag, d.P / psi, string(d.route)}; %#ok<SAGROW>
    fprintf('%-9s Pdew = %.2f psia\n', s.tag, d.P / psi);

    P = Pall;
    n = numel(P);
    [ph, SL, rL, rV, xC8, yC8, ZL, ZV] = deal(nan(n, 1));
    ok = false(n, 1);
    K = d.K;                                        % dew-point K seeds the first two-phase flash
    for i = 1:n
        r = pt.solve(P(i) * psi, cfg.T, cfg.z, s.r, 'K_seed', K);
        ok(i) = r.converged;
        if ~r.converged, continue; end
        ph(i) = r.phases;
        if r.phases == 1
            SL(i) = 0;                              % dew-point fluid above Pdew: gas only
        else
            SL(i) = r.SL;  K = r.K;
        end
        rL(i) = r.rhoM_L;  rV(i) = r.rhoM_V;         % kg/m3 (PTFlash: 1e-3*MW/V)
        xC8(i) = r.x(2);   yC8(i) = r.y(2);
        ZL(i) = r.Z_L;     ZV(i) = r.Z_V;
    end
    out = table(P, ph, SL, rL, rV, xC8, yC8, ZL, ZV, ok, 'VariableNames', ...
        {'P_psia', 'phases', 'SL', 'rhoL_kgm3', 'rhoV_kgm3', 'x_nC8', 'y_nC8', 'Z_L', 'Z_V', 'converged'});
    out = sortrows(out, 'P_psia');
    writetable(out, fullfile(cfg.dir, "matlab_" + s.tag + ".csv"));
    fprintf('          %d pressures, %d not converged\n', n, sum(~ok));
end
writetable(dew, fullfile(cfg.dir, "matlab_dewpoints.csv"));
