% =========================================================================
% DYNAMIC CO2 CONFINED TIP BACK-CALIBRATION UTILITY
% Inverts confined dew-point experimental targets directly from Mix_Props
% to extract exact kijc values and energy disparity ratios (E*)
% =========================================================================
clc; clear; close all;
format shortG;

% Add src path if running from tools/ folder or root
if isfolder('../src'), addpath('../src'); end
if isfolder('src'), addpath('src'); end

xlsxPath = 'config/MixtureData.xlsx';
if ~isfile(xlsxPath)
    xlsxPath = '../config/MixtureData.xlsx';
    if ~isfile(xlsxPath)
        xlsxPath = 'MixtureData_2.xlsx';
        if ~isfile(xlsxPath), xlsxPath = '../MixtureData_2.xlsx'; end
    end
end

fprintf('===================================================================\n');
fprintf('     AUTOMATED CO2 CONFINED TIP BACK-CALIBRATION ENGINE            \n');
fprintf('===================================================================\n\n');

%% 1. READ AND FILTER EXPERIMENTAL DATASETS FROM EXCEL
optsMix = detectImportOptions(xlsxPath, 'Sheet', 'Mix_Props');
optsMix.VariableNamingRule = 'preserve';
mixTable = readtable(xlsxPath, optsMix);

% Find all rows containing CO2 in the Mixture column
co2Mask = contains(string(mixTable.("Mixture")), "CO2");
co2Rows = mixTable(co2Mask, :);

if isempty(co2Rows)
    error('DataError:NoCO2Data', 'No CO2-included mixture rows found in the Mix_Props sheet of %s', xlsxPath);
end

fprintf('Available CO2 Experimental Datasets in database:\n');
nCO2 = height(co2Rows);
for idx = 1:nCO2
    fprintf('  %d. Mixture: %-12s | Rock: %-4s | T: %5.1f°C | r_p: %4.1f nm | Exp Confined Pdew: %g psi\n', ...
        idx, string(co2Rows.("Mixture")(idx)), string(co2Rows.("Rock")(idx)), ...
        co2Rows.("Temp_C")(idx), co2Rows.("Pore_Radius_nm")(idx), co2Rows.("Confined_Pdew")(idx));
end
fprintf('  %d. [ALL DATASETS] -> Run batch calibration report across all CO2 rows\n', nCO2 + 1);

choice = input(sprintf('\nSelect dataset index to calibrate (1-%d): ', nCO2 + 1));
if isempty(choice) || choice < 1 || choice > (nCO2 + 1), choice = 1; end

if choice <= nCO2
    runIndices = choice;
else
    runIndices = 1:nCO2;
end

%% 2. EXECUTE ROOT INVERSION LOOP
psi2Pa = 6894.757;
resultsSummary = struct();

for runIdx = 1:length(runIndices)
    rowIdx   = runIndices(runIdx);
    mixName  = char(string(co2Rows.("Mixture")(rowIdx)));
    rockName = char(string(co2Rows.("Rock")(rowIdx)));
    T_C      = co2Rows.("Temp_C")(rowIdx);
    T_K      = T_C + 273.15;
    r_nm     = co2Rows.("Pore_Radius_nm")(rowIdx);
    r_cap    = r_nm * 1e-9;
    
    P_exp_psi = co2Rows.("Confined_Pdew")(rowIdx);
    P_exp_Pa  = P_exp_psi * psi2Pa;
    
    % Parse feed composition string e.g. "[0.100, 0.828, 0.072]"
    rawStr   = string(co2Rows.("mole_frac")(rowIdx));
    cleanStr = strrep(strrep(rawStr, "[", ""), "]", "");
    z_feed   = str2double(strsplit(cleanStr, ","));
    z_feed   = z_feed(:) / sum(z_feed);
    
    fprintf('\n-------------------------------------------------------------------\n');
    fprintf(' Running Calibration [%d/%d]: Mixture = %s | Rock = %s\n', ...
        runIdx, length(runIndices), mixName, rockName);
    fprintf(' Conditions: T = %.2f°C (%.2f K) | r_p = %.1f nm | Target Pdew = %.2f psi\n', ...
        T_C, T_K, r_nm, P_exp_psi);
    fprintf(' Feed Composition: z = [%s]\n', num2str(z_feed'));
    
    % Load Entities
    fluid = entities.FluidProperties.loadFromWorkbook(xlsxPath, mixName);
    rock  = entities.RockProperties.loadFromWorkbook(xlsxPath, rockName);
    nc = fluid.NC; NA = 6.02214076e23;
    
    % Locate CO2 and non-C1 Heavy Alkane indices
    idx_CO2    = find(strcmpi(fluid.ComponentNames, 'CO2'), 1);
    idx_Alkane = find(~strcmpi(fluid.ComponentNames, 'CO2') & ~strcmpi(fluid.ComponentNames, 'C1'), 1);
    
    if isempty(idx_CO2) || isempty(idx_Alkane)
        error('ComponentError:MappingFailed', 'Could not identify CO2 or Heavy Alkane component in [%s]', mixName);
    end
    
    alkaneName = fluid.ComponentNames{idx_Alkane};
    fprintf(' Calibrating Interaction Pair: CO2 (%d) <-> %s (%d)\n', idx_CO2, alkaneName, idx_Alkane);
    
    % Energy calculations for E* ratio
    eps_fluid_molar = fluid.LJ_Energy(:) * NA;
    eps_rock_molar  = rock.EnergyVector(:);
    w_norm          = rock.WeightVector(:);
    
    eps_wall_comp = zeros(nc, 1);
    for i = 1:nc
        eps_wall_comp(i) = sum(w_norm .* sqrt(eps_rock_molar .* eps_fluid_molar(i)));
    end
    
    E_star_val = sqrt(eps_fluid_molar(idx_CO2) * eps_fluid_molar(idx_Alkane)) / ...
                 sqrt(eps_wall_comp(idx_CO2) * eps_wall_comp(idx_Alkane));
             
    % Setup Root-Finding Residual Function
    targetResidual = @(k_val) evaluatePdewError(...
        k_val, idx_CO2, idx_Alkane, fluid, rock, T_K, z_feed, r_cap, P_exp_Pa);
    
    % Execute fzero search bounded between -0.8 and 0.8
    options = optimset('TolX', 1e-5, 'Display', 'off');
    try
        k_co2_opt = fzero(targetResidual, [-0.8, 0.8], options);
        
        % Store Results
        resultsSummary(runIdx).Mixture      = mixName;
        resultsSummary(runIdx).Rock         = rockName;
        resultsSummary(runIdx).Pair         = sprintf('CO2-%s', alkaneName);
        resultsSummary(runIdx).Temp_C       = T_C;
        resultsSummary(runIdx).Radius_nm    = r_nm;
        resultsSummary(runIdx).Target_Pdew  = P_exp_psi;
        resultsSummary(runIdx).Calibrated_K = k_co2_opt;
        resultsSummary(runIdx).E_star       = E_star_val;
        
        fprintf('  [SUCCESS] Calibrated kijc (CO2-%s) = %10.6f | E* = %10.6f\n', ...
            alkaneName, k_co2_opt, E_star_val);
    catch ME
        fprintf('  [FAILED] Inversion did not converge: %s\n', ME.message);
    end
end

%% 3. DISPLAY FINAL CALIBRATION REPORT
if ~isempty(fieldnames(resultsSummary))
    fprintf('\n===================================================================\n');
    fprintf('                  CO2 TIP CALIBRATION SUMMARY TABLE                \n');
    fprintf('===================================================================\n');
    fprintf(' %-12s | %-4s | %-10s | %-6s | %-6s | %-10s | %-10s\n', ...
        'Mixture', 'Rock', 'Pair', 'T [°C]', 'r [nm]', 'Calib kijc', 'Ratio E*');
    fprintf('-------------------------------------------------------------------\n');
    for i = 1:length(resultsSummary)
        fprintf(' %-12s | %-4s | %-10s | %6.1f | %6.1f | %10.6f | %10.6f\n', ...
            resultsSummary(i).Mixture, resultsSummary(i).Rock, resultsSummary(i).Pair, ...
            resultsSummary(i).Temp_C, resultsSummary(i).Radius_nm, ...
            resultsSummary(i).Calibrated_K, resultsSummary(i).E_star);
    end
    fprintf('===================================================================\n');
end

% --- Helper Function for Root Finding ---
function err = evaluatePdewError(k_co2_alkane, idx1, idx2, fluid, rock, T, z, r_cap, P_target)
    nc = fluid.NC;
    kijc = zeros(nc, nc);
    kijc(idx1, idx2) = k_co2_alkane; 
    kijc(idx2, idx1) = k_co2_alkane;
    
    eos = thermo.ConfinedEOS(fluid, rock, 'kijc', kijc);
    stability = solvers.StabilityTester(eos);
    flash = solvers.FlashEngine(eos, stability);
    
    [P_pred, ~, ~, ~, ~] = flash.solveDewPoint(T, P_target, z, r_cap, 'Solver', 'newton');
    err = P_pred - P_target;
end