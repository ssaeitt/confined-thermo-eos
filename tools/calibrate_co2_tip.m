% =========================================================================
% DYNAMIC CO2 CONFINED TIP BACK-CALIBRATION UTILITY (SIMPLEX OPTIMIZER)
% Uses Nelder-Mead Simplex (fminsearch) to invert confined dew-point 
% experimental targets and extract exact negative TIP values.
% =========================================================================
clc; clear; close all;
format shortG;

% Add src path if running from tools/ folder or root
if isfolder('../src'), addpath('../src'); end
if isfolder('src'), addpath('src'); end

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
        'Could not locate MixtureData.xlsx in config/ or root directory.');
end

fprintf('===================================================================\n');
fprintf('     AUTOMATED CO2 CONFINED TIP BACK-CALIBRATION ENGINE            \n');
fprintf('===================================================================\n\n');
fprintf(' Using Database File: %s\n\n', xlsxPath);

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

% Extract available column names
varNames = string(co2Rows.Properties.VariableNames);
nCO2 = height(co2Rows);

% Pre-extract row parameters with robust fallback logic
parsedData = struct();
for idx = 1:nCO2
    mixName  = char(string(co2Rows.("Mixture")(idx)));
    rockName = char(string(co2Rows.("Rock")(idx)));
    
    % 1. Target Confined Dew Point [psi]
    if ismember("Confined_Pdew", varNames)
        P_exp = co2Rows.("Confined_Pdew")(idx);
    elseif ismember("Pdew_Exp_upper", varNames)
        P_exp = co2Rows.("Pdew_Exp_upper")(idx);
    else
        P_exp = co2Rows.("Bulk_Pdew")(idx);
    end
    
    % 2. Temperature [°C] Fallback Mapping
    if ismember("Temp_C", varNames) && ~isnan(co2Rows.("Temp_C")(idx))
        T_C = co2Rows.("Temp_C")(idx);
    elseif contains(mixName, "nC8")
        T_C = 20.0;
    elseif contains(mixName, "nC5")
        T_C = 37.8;
    else
        T_C = 20.0;
    end
    
    % 3. Pore Radius [nm] Fallback Mapping
    if ismember("Pore_Radius_nm", varNames) && ~isnan(co2Rows.("Pore_Radius_nm")(idx))
        r_nm = co2Rows.("Pore_Radius_nm")(idx);
    elseif strcmpi(rockName, "EF2") || contains(mixName, "nC8")
        r_nm = 54.0;  % EF2 matrix average pore radius
    elseif strcmpi(rockName, "B1") || contains(mixName, "nC5")
        r_nm = 11.25; % B1 nanopore radius
    else
        r_nm = 10.0;
    end
    
    parsedData(idx).Mixture     = mixName;
    parsedData(idx).Rock        = rockName;
    parsedData(idx).Temp_C      = T_C;
    parsedData(idx).Radius_nm   = r_nm;
    parsedData(idx).Target_Pdew = P_exp;
    parsedData(idx).MoleFrac    = string(co2Rows.("mole_frac")(idx));
end

%% 2. EXECUTE ROOT INVERSION LOOP
psi2Pa = 6894.757;
resultsSummary = struct();

for runIdx = 1:nCO2
    mixName  = parsedData(runIdx).Mixture;
    rockName = parsedData(runIdx).Rock;
    T_C      = parsedData(runIdx).Temp_C;
    T_K      = T_C + 273.15;
    r_nm     = parsedData(runIdx).Radius_nm;
    r_cap    = r_nm * 1e-9;
    
    P_exp_psi = parsedData(runIdx).Target_Pdew;
    P_exp_Pa  = P_exp_psi * psi2Pa;
    
    % Parse feed composition string e.g. "[0.100, 0.828, 0.072]"
    rawStr   = parsedData(runIdx).MoleFrac;
    cleanStr = strrep(strrep(rawStr, "[", ""), "]", "");
    z_feed   = str2double(strsplit(cleanStr, ","));
    z_feed   = z_feed(:) / sum(z_feed);
    
    fprintf('\n-------------------------------------------------------------------\n');
    fprintf(' Running Calibration [%d/%d]: Mixture = %s | Rock = %s\n', ...
        runIdx, nCO2, mixName, rockName);
    fprintf(' Conditions: T = %.2f°C (%.2f K) | r_p = %.1f nm | Target Pdew = %.2f psi\n', ...
        T_C, r_nm, P_exp_psi);
    
    % Load Entities
    fluid = entities.FluidProperties.loadFromWorkbook(xlsxPath, mixName);
    rock  = entities.RockProperties.loadFromWorkbook(xlsxPath, rockName);
    nc = fluid.NC; NA = 6.02214076e23;
    
    % Locate CO2 and non-C1 Heavy Alkane indices
    idx_CO2    = find(strcmpi(fluid.ComponentNames, 'CO2'), 1);
    idx_C1     = find(strcmpi(fluid.ComponentNames, 'C1'), 1);
    idx_Alkane = find(~strcmpi(fluid.ComponentNames, 'CO2') & ~strcmpi(fluid.ComponentNames, 'C1'), 1);
    
    alkaneName = fluid.ComponentNames{idx_Alkane};
    
    % Energy calculations for E* ratio
    eps_fluid_molar = fluid.LJ_Energy(:) * NA;
    eps_rock_molar  = rock.EnergyVector(:);
    w_norm          = rock.WeightVector(:);
    
    eps_wall_comp = zeros(nc, 1);
     for i = 1:nc
        eps_wall_comp(i) = sum(w_norm .* sqrt(eps_rock_molar .* eps_fluid_molar(i)));
    end
    
    E_star_co2_alkane = sqrt(eps_fluid_molar(idx_CO2) * eps_fluid_molar(idx_Alkane)) / ...
                        sqrt(eps_wall_comp(idx_CO2) * eps_wall_comp(idx_Alkane));
                    
    E_star_c1_alkane = sqrt(eps_fluid_molar(idx_C1) * eps_fluid_molar(idx_Alkane)) / ...
                       sqrt(eps_wall_comp(idx_C1) * eps_wall_comp(idx_Alkane));
             
    % Negative Hydrocarbon TIP baseline (e.g. -1.0 to -2.0 depending on E*)
    k_c1_alkane_c = -0.35 * E_star_c1_alkane;

    fprintf(' Hydrocarbon TIP Baseline (C1-%s): TIP = %.5f\n', alkaneName, k_c1_alkane_c);
    fprintf(' Optimizing CO2 TIP Pair via Simplex: CO2 <-> %s\n', alkaneName);

    % Objective Function for Nelder-Mead Simplex (Absolute Pressure Residual in psi)
    objFun = @(k_val) evaluateObjective(...
        k_val, idx_CO2, idx_Alkane, idx_C1, k_c1_alkane_c, ...
        fluid, rock, T_K, z_feed, r_cap, P_exp_Pa);
    
    % Set initial guess based on mixture type
    if contains(mixName, "nC8")
        k_init = -4.0;
    else
        k_init = -0.8;
    end
    
    options = optimset('TolX', 1e-4, 'TolFun', 1e-2, 'Display', 'iter', 'MaxIter', 100);
    
    try
        [k_co2_opt, fval] = fminsearch(objFun, k_init, options);
        
        % Calculate final predicted pressure
        eos_final = thermo.ConfinedEOS(fluid, rock, 'kijc', buildKijc(k_co2_opt, idx_CO2, idx_Alkane, idx_C1, k_c1_alkane_c, nc));
        flash_final = solvers.FlashEngine(eos_final, solvers.StabilityTester(eos_final));
        [P_pred_Pa, ~, ~, ~, ~] = flash_final.solveDewPoint(T_K, P_exp_Pa, z_feed, r_cap, 'Solver', 'newton');
        P_pred_psi = P_pred_Pa / psi2Pa;
        
        resultsSummary(runIdx).Mixture      = mixName;
        resultsSummary(runIdx).Rock         = rockName;
        resultsSummary(runIdx).Pair         = sprintf('CO2-%s', alkaneName);
        resultsSummary(runIdx).Temp_C       = T_C;
        resultsSummary(runIdx).Radius_nm    = r_nm;
        resultsSummary(runIdx).Target_Pdew  = P_exp_psi;
        resultsSummary(runIdx).Pred_Pdew    = P_pred_psi;
        resultsSummary(runIdx).Calibrated_K = k_co2_opt;
        resultsSummary(runIdx).E_star       = E_star_co2_alkane;
        
        fprintf('\n  [SUCCESS] Calibrated TIP (CO2-%s) = %10.6f | Pred Pdew = %.2f psi | Target = %.2f psi | E* = %10.6f\n\n', ...
            alkaneName, k_co2_opt, P_pred_psi, P_exp_psi, E_star_co2_alkane);
    catch ME
        fprintf('\n  [FAILED] Simplex optimization did not converge: %s\n\n', ME.message);
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
            resultsSummary(i).Target_Pdew, resultsSummary(i).Pred_Pdew, resultsSummary(i).Calibrated_K);
    end
    fprintf('===================================================================\n');
end

% --- Helper Function for Root Finding ---
function absErr_psi = evaluateObjective(k_co2_alkane, idx_CO2, idx_Alkane, idx_C1, k_c1_alkane_c, fluid, rock, T, z, r_cap, P_target)
    nc = fluid.NC;
    kijc = buildKijc(k_co2_alkane, idx_CO2, idx_Alkane, idx_C1, k_c1_alkane_c, nc);
    
    eos = thermo.ConfinedEOS(fluid, rock, 'kijc', kijc);
    stability = solvers.StabilityTester(eos);
    flash = solvers.FlashEngine(eos, stability);
    
    [P_pred, ~, ~, ~, stats] = flash.solveDewPoint(T, P_target, z, r_cap, 'Solver', 'newton');
    if ~stats.converged || isnan(P_pred)
        absErr_psi = 1e6; % Penalty for non-convergence
    else
        absErr_psi = abs(P_pred - P_target) / 6894.757; % Minimize absolute error in psi
    end
end

function kijc = buildKijc(k_co2_alkane, idx_CO2, idx_Alkane, idx_C1, k_c1_alkane_c, nc)
    kijc = zeros(nc, nc);
    kijc(idx_CO2, idx_Alkane) = k_co2_alkane;
    kijc(idx_Alkane, idx_CO2) = k_co2_alkane;
    
    if ~isempty(idx_C1) && ~isempty(idx_Alkane)
        kijc(idx_C1, idx_Alkane) = k_c1_alkane_c;
        kijc(idx_Alkane, idx_C1) = k_c1_alkane_c;
    end
end