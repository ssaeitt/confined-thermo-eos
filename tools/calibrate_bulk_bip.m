% =========================================================================
% AUTOMATED BULK BIP (k_ij) OPTIMIZER FOR CO2-ALKANE PAIRS
% Solves for exact bulk k_ij values matching experimental bulk Pdew
% =========================================================================
clc; clear; close all;

% Smart path resolution for src directory
if isfolder('../src'), addpath('../src'); end
if isfolder('src'), addpath('src'); end

% 2. Smart Path Resolution for Database Excel File
possiblePaths = { ...
    'config/MixtureData.xlsx', ...
    '../config/MixtureData.xlsx', ...
    'MixtureData_2.xlsx', ...
    '../MixtureData_2.xlsx' ...
};

xlsxPath = '';
for i = 1:length(possiblePaths)
    if isfile(possiblePaths{i})
        xlsxPath = possiblePaths{i};
        break;
    end
end

if isempty(xlsxPath)
    error('FileIOError:SpreadsheetNotFound', ...
        'Could not locate MixtureData.xlsx or MixtureData_2.xlsx in root or config/ folders.');
end

psi2Pa = 6894.757;

% Define Bulk Calibration Targets
targets = struct(...
    'Mix',    {'CO2-C1-nC8', 'CO2-C1-nC5'}, ...
    'Rock',   {'EF2', 'B1'}, ...
    'Temp_C', {20.0, 37.8}, ...
    'z',      {[0.100; 0.828; 0.072], [0.150; 0.700; 0.150]}, ...
    'P_exp',  {4126, 2701});

fprintf('===================================================================\n');
fprintf('         BULK BIP (k_ij) INVERSION RUNNER FOR CO2 SYSTEMS          \n');
fprintf('===================================================================\n\n');
fprintf(' Using Database File: %s\n\n', xlsxPath);

for i = 1:length(targets)
    mixName = targets(i).Mix;
    T_K     = targets(i).Temp_C + 273.15;
    z_feed  = targets(i).z;
    P_exp_Pa = targets(i).P_exp * psi2Pa;
    
    fluid = entities.FluidProperties.loadFromWorkbook(xlsxPath, mixName);
    rock  = entities.RockProperties.loadFromWorkbook(xlsxPath, targets(i).Rock);
    
    idx_CO2    = find(strcmpi(fluid.ComponentNames, 'CO2'), 1);
    idx_Alkane = find(~strcmpi(fluid.ComponentNames, 'CO2') & ~strcmpi(fluid.ComponentNames, 'C1'), 1);
    alkaneName = fluid.ComponentNames{idx_Alkane};
    
    fprintf('Calibrating Pair: CO2-%-4s (Target = %4d psi)...\n', alkaneName, targets(i).P_exp);

    % Residual function tuning BIP at r_cap = Inf (Bulk)
    resFun = @(k_val) runBulkError(k_val, idx_CO2, idx_Alkane, fluid, rock, T_K, z_feed, P_exp_Pa);
    
    options = optimset('TolX', 1e-4, 'Display', 'off');

    % Test bracket endpoints explicitly before fzero
    k_min = 0.05; k_max = 0.50;
    err_min = resFun(k_min);
    err_max = resFun(k_max);
    
    if (err_min * err_max) > 0
        % If same sign, try expanding search bounds
        k_min = -0.10; k_max = 0.60;
        err_min = resFun(k_min);
        err_max = resFun(k_max);
    end
    
    try
        k_opt = fzero(resFun, [k_min, k_max], options);
        fprintf('  [SUCCESS] Mixture: %-12s | Pair: CO2-%-4s | Target: %4d psi | Calibrated Bulk k_ij: %.4f\n\n', ...
            mixName, alkaneName, targets(i).P_exp, k_opt);
    catch ME
        fprintf('  [FAILED] Mixture: %-12s | Pair: CO2-%-4s | Target: %4d psi | Bounds [%.2f, %.2f] | Error: %s\n\n', ...
            mixName, alkaneName, targets(i).P_exp, k_min, k_max, ME.message);
    end
end

function err = runBulkError(k_val, idx1, idx2, fluid, rock, T, z, P_target)
    % Extract existing BIPMatrix and update the target pair symmetrically
    bip = fluid.BIPMatrix;
    bip(idx1, idx2) = k_val;
    bip(idx2, idx1) = k_val;
    
    % Re-instantiate FluidProperties with the updated BIPMatrix
    fluid_mod = entities.FluidProperties(...
        fluid.ComponentNames, fluid.Tc, fluid.Pc, fluid.omega, fluid.MW, ...
        fluid.Parachor, fluid.VolumeShift, fluid.LJ_Size, fluid.LJ_Energy, bip);
    
    % Pass updated fluid entity to ConfinedEOS
    eos = thermo.ConfinedEOS(fluid_mod, rock, 'kijc', zeros(fluid_mod.NC, fluid_mod.NC));
    stability = solvers.StabilityTester(eos);
    flash = solvers.FlashEngine(eos, stability);
    
    % Unconfined Bulk Flash (r_cap = Inf)
    [P_pred, ~, ~, ~, stats] = flash.solveDewPoint(T, P_target, z, Inf, 'Solver', 'newton');

    if ~stats.converged || isnan(P_pred)
        err = 1e8; % Penalize non-convergence
    else
        err = P_pred - P_target;
    end
end