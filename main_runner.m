% =========================================================================
% CONFINED THERMODYNAMIC EOS MODELING: INTERACTIVE PDEW DRIVER
% Object-Oriented framework for VLE under nanoporous confinement
% Supports both Combined (FWI + Pcap) and FWI-Only (Pcap = 0) physics modes
% =========================================================================
clc; clear; close all;
format shortG;

% 1. INITIALIZE WORKSPACE & NAMESPACE PATHS
if ~isfolder('src')
    error('ArchitectureError:MissingSource', ...
        'The "src" directory was not found. Ensure you are running from the repository root.');
end
addpath('src');

possiblePaths = { ...
    'config/MixtureData.xlsx', ...
    '../config/MixtureData.xlsx', ...
    'MixtureData.xlsx', ...
    '../MixtureData.xlsx', ...
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
    error('DataIOError:WorkbookNotFound', ...
        'Thermodynamic database spreadsheet not found in config/ or root directory.');
end

fprintf('===================================================================\n');
fprintf('   NANOPOROUS CONFINEMENT PDEW PREDICTION ENGINE (OOP V2.3)        \n');
fprintf('===================================================================\n\n');

%% 2. INTERACTIVE RESEARCH PARAMETER & MIXTURE CONFIGURATION
% --- Model Confinement Physics Selection ---
% --- Model Confinement Physics Selection ---
fprintf('Select Confinement Physics Mode:\n');
fprintf('  1. Combined Model (Fluid-Wall Interaction + Pcap)\n');
fprintf('  2. FWI-Only Model (Fluid-Wall Interaction Only, Pcap = 0)\n');
modeChoice = input('Enter physics mode index (1-2, default = 2): ');
if isempty(modeChoice) || ~ismember(modeChoice, [1, 2]), modeChoice = 2; end
if modeChoice == 1
    modelMode = 'Combined';
else
    modelMode = 'FWI-Only';
end
fprintf('Selected Physics Mode: [%s]\n\n', modelMode);

% --- Numerical Solver Algorithm Selection ---
fprintf('Select Numerical Solver Engine:\n');
fprintf('  1. Newton-Raphson (Full Jacobian + Armijo Line Search)\n');
fprintf('  2. Quasi-Newton Broyden (Rank-1 Sherman-Morrison Update)\n');
fprintf('  3. Successive Substitution (SS Standalone Fixed-Point Solver)\n');
solverChoice = input('Enter solver index (1-3, default = 1): ');
if isempty(solverChoice) || ~ismember(solverChoice, 1:3), solverChoice = 1; end

solverMap = {'newton', 'quasinewton', 'ss'};
solverNameMap = {'Newton-Raphson', 'Quasi-Newton (Broyden)', 'Successive Substitution (SS)'};
selectedSolver = solverMap{solverChoice};
selectedSolverName = solverNameMap{solverChoice};

% --- Successive Substitution (SS) Preconditioner Control ---
if strcmpi(selectedSolver, 'ss')
    useSSPrecond = true;
    fprintf('Selected Main Solver: [%s]\n\n', selectedSolverName);
else
    fprintf('Enable Successive Substitution (SS) Preconditioning?\n');
    fprintf('  1. Yes (Recommended: Runs 30 SS coarse steps before Newton)\n');
    fprintf('  2. No  (Direct gradient search from Wilson/TPD seed)\n');
    precondChoice = input('Enter preconditioning choice (1-2, default = 1): ');
    if isempty(precondChoice) || ~ismember(precondChoice, [1, 2]), precondChoice = 1; end
    useSSPrecond = (precondChoice == 1);
    fprintf('Selected Main Solver: [%s] | SS Preconditioner: [%s]\n\n', ...
        selectedSolverName, mat2str(useSSPrecond));
end

% --- Geological Formation Selection ---
fprintf('Select Target Geological Formation:\n');
fprintf('  1. Eagle Ford (EF2  - Carbonate-Rich, Low TOC)\n');
fprintf('  2. Barnett    (B1   - Silicate-Rich,  High TOC)\n');
rockChoice = input('Enter formation index (1-2, default = 2): ');
rockMap = {'EF2', 'B1'};
if isempty(rockChoice) || ~ismember(rockChoice, [1, 2]), rockChoice = 2; end
selectedRock = rockMap{rockChoice};

% --- Hydrocarbon & CO2 Mixture System Selection ---
fprintf('\nSelect Gas Mixture System:\n');
fprintf('  1. C1-nC5      (Methane / n-Pentane)\n');
fprintf('  2. C1-nC8      (Methane / n-Octane)\n');
fprintf('  3. C1-nC10     (Methane / n-Decane)\n');
fprintf('  4. C1-nC5-nC10 (Ternary Hydrocarbon System)\n');
fprintf('  5. CO2-C1-nC8  (Carbon Dioxide / Methane / n-Octane)\n');
fprintf('  6. CO2-C1-nC5  (Carbon Dioxide / Methane / n-Pentane)\n');
mixChoice = input('Enter mixture index (1-6, default = 2): ');
mixMap = {'C1-nC5', 'C1-nC8', 'C1-nC10', 'C1-nC5-nC10', 'CO2-C1-nC8', 'CO2-C1-nC5'};
if isempty(mixChoice) || ~ismember(mixChoice, 1:6), mixChoice = 2; end
selectedMixName = mixMap{mixChoice};

% --- Thermodynamic State & Confinement Boundaries ---
defaultT_C = 67.85;
defaultP_MPa = 25.0;

if strcmpi(selectedMixName, 'C1-nC8')
    defaultT_C = 20.0;
    defaultP_MPa = 28.0; % Baseline near 4000 psi
elseif strcmpi(selectedMixName, 'CO2-C1-nC8')
    defaultT_C = 20.0;
    defaultP_MPa = 23.0; % Baseline near 3300 psi
elseif strcmpi(selectedMixName, 'CO2-C1-nC5')
    defaultT_C = 37.8;
    defaultP_MPa = 15.0;
end

T_C = input(sprintf('\nEnter system temperature [°C] (default = %.2f): ', defaultT_C));
if isempty(T_C), T_C = defaultT_C; end
T_K = T_C + 273.15;

default_r_nm = 11.25;
if strcmpi(selectedRock, 'EF2'), default_r_nm = 54.0; end

r_nm = input(sprintf('Enter pore radius [nm] (default = %.2f, enter Inf for bulk): ', default_r_nm));
if isempty(r_nm), r_nm = default_r_nm; end
if isinf(r_nm)
    r_cap = Inf;
else
    r_cap = r_nm * 1e-9;
end

P_guess_MPa = input(sprintf('Enter initial Pdew guess [MPa] (default = %.1f): ', defaultP_MPa));
if isempty(P_guess_MPa), P_guess_MPa = defaultP_MPa; end
P_guess_Pa = P_guess_MPa * 1e6;

%% 3. LOAD DATA MODELS & RESOLVE MIXTURE COMPOSITIONS
fprintf('\n[1/4] Ingesting database entities from workbook: %s...\n', xlsxPath);

% Load fluid properties filtered specifically for the chosen mixture components
fluid = entities.FluidProperties.loadFromWorkbook(xlsxPath, selectedMixName);
rock  = entities.RockProperties.loadFromWorkbook(xlsxPath, selectedRock);

% --- Interactive Composition Selector from Mix_Props ---
optsMix = detectImportOptions(xlsxPath, 'Sheet', 'Mix_Props');
optsMix.VariableNamingRule = 'preserve';
mixTable = readtable(xlsxPath, optsMix);

% Filter table for active mixture and rock formation
validRows = mixTable(strcmp(mixTable.("Mixture"), selectedMixName) & ...
                     strcmp(mixTable.("Rock"), selectedRock), :);
fprintf('\nAvailable Experimental Datasets for [%s] in [%s]:\n', selectedMixName, selectedRock);
nRows = height(validRows);
for idx = 1:nRows
    fprintf('  %d. Feed z = %s | Exp Bulk Pdew = %g psi\n', ...
        idx, string(validRows.("mole_frac")(idx)), validRows.("Bulk_Pdew")(idx));
end
fprintf('  %d. [Custom Input] -> Enter custom molar fractions manually\n', nRows + 1);

compChoice = input(sprintf('Select feed composition index (1-%d, default = 1): ', nRows + 1));
if isempty(compChoice) || compChoice < 1 || compChoice > (nRows + 1)
    compChoice = 1;
end

hasExpMatch = false;
expRow = [];
if compChoice <= nRows
    hasExpMatch = true;
    expRow = validRows(compChoice, :);
    rawStr = string(validRows.("mole_frac")(compChoice));
    cleanStr = strrep(strrep(rawStr, "[", ""), "]", "");
    z_feed = str2double(strsplit(cleanStr, ","));
    z_feed = z_feed(:) / sum(z_feed);
    fprintf('      -> Selected Dataset Feed: z = [%s]\n', num2str(z_feed', '%.4f '));
else
    % Prompt user for custom vector
    fprintf('Enter molar composition matching [%s] (e.g., [0.90, 0.10]):\n', ...
        strjoin(fluid.ComponentNames, ', '));
    z_custom = input('z = ');
    if isempty(z_custom) || length(z_custom) ~= fluid.NC
        error('InputError:InvalidComposition', 'Vector must contain exactly %d elements.', fluid.NC);
    end
    z_feed = z_custom(:) / sum(z_custom);
    fprintf('      -> Applied Custom Feed: z = [%s]\n', num2str(z_feed', '%.4f '));
end

z_feed = z_feed(:);

%% 4. EXPLICIT TIP MATRIX (kijc) VIA INVERSE-DISPARITY CORRELATION
fprintf('[2/4] Computing Ternary Interaction Parameter (TIP) matrix...\n');

% Common EOS FWI Tuners
k_val       = 64.4376;
pT_wall_val = 1.3745;
lambda_val  = 10.0338;

% 1. Define Correlation Fitting Coefficients (Inverse-Disparity Model)
if strcmpi(modelMode, 'FWI-Only')
    A_corr =  101.415477921051;
    B_corr =    1.23381011476174;
    C_corr = -168.946183223101;
else
    A_corr =   90.5402886251168;
    B_corr =    1.09588610509636;
    C_corr = -150.825648642323;
end

nc = fluid.NC;
NA = 6.02214076e23; % Avogadro's Constant [1/mol]

% 2. Extract and Convert Energy Profiles to Molar Basis [J/mol]
eps_fluid_molar = fluid.LJ_Energy(:) * NA;
eps_rock_molar  = rock.EnergyVector(:);
w_norm          = rock.WeightVector(:);
rawTOC          = max(rock.RawTOC, 1e-4);

% 3. Compute Composite Molar Fluid-Wall Energy per Component
eps_wall_comp = zeros(nc, 1);
for i = 1:nc
    % Geometric mean mixing between fluid component i and each mineral group
    eps_wall_comp(i) = sum(w_norm .* sqrt(eps_rock_molar .* eps_fluid_molar(i)));
end

% 4. Populate Symmetric TIP Matrix (kijc)
kijc_matrix = zeros(nc, nc);
for r = 1:nc
    for c = (r+1):nc
        % Dimensionless energy disparity ratio
        E_star = sqrt(eps_fluid_molar(r) * eps_fluid_molar(c)) / ...
                 sqrt(eps_wall_comp(r) * eps_wall_comp(c));
        
        % Inverse-Disparity empirical correlation
        val_corr = A_corr * E_star + (B_corr * E_star) / rawTOC + C_corr;
        
        kijc_matrix(r, c) = val_corr;
        kijc_matrix(c, r) = val_corr;
    end
end

fprintf('      -> TIP Matrix Configured (Inverse-Disparity Model | TOC = %.2f%%)\n', rock.w_TOC * 100);

%% 5. INSTANTIATE NUMERICAL SOLVER ENGINES
fprintf('[3/4] Initializing Unconfined Bulk Engine (kijc = 0)...\n');
eos_bulk       = thermo.ConfinedEOS(fluid, rock, 'kijc', zeros(nc, nc));
stability_bulk = solvers.StabilityTester(eos_bulk, 'MaxIterations', 500, 'GradTolerance', 1e-6);
flash_bulk     = solvers.FlashEngine(eos_bulk, stability_bulk, ...
    'MaxIterations', 1000, ...
    'Tolerance', 1e-6, ...
    'PreconditionWithSS', useSSPrecond);

fprintf('[4/4] Initializing Confined Confinement Engine (Active kijc)...\n');
eos_conf       = thermo.ConfinedEOS(fluid, rock, 'kijc', kijc_matrix);
stability_conf = solvers.StabilityTester(eos_conf, 'MaxIterations', 500, 'GradTolerance', 1e-6);
flash_conf     = solvers.FlashEngine(eos_conf, stability_conf, ...
    'MaxIterations', 1000, ...
    'Tolerance', 1e-6, ...
    'PreconditionWithSS', useSSPrecond);

%% 6. EXECUTE TWO-STAGE SATURATION POINT MOLECULAR CONTINUATION
tic;

% --- STAGE 1: UNCONFINED BULK CALIBRATION ANCHOR ---
fprintf('\nExecuting Stage 1: Solving Unconfined Bulk Saturation Boundary...\n');
[P_bulk_Pa, K_bulk, ~, ~, bulk_stats] = flash_bulk.solveDewPoint(...
    T_K, P_guess_Pa, z_feed, Inf, ...
    'Solver', selectedSolver, ...
    'CapMode', 'PL', ...
    'PreconditionSS', useSSPrecond);

psi2Pa = 6894.757;
if bulk_stats.converged
    fprintf('      -> Bulk Stage Converged: %.4f MPa (%8.2f psia) in %d iter (Norm = %.2e)\n', ...
        P_bulk_Pa / 1e6, P_bulk_Pa / psi2Pa, bulk_stats.iterations, bulk_stats.residual_norm);
else
    fprintf('      -> Bulk Stage Stalled (%d iter, norm = %.2e).\n', ...
        bulk_stats.iterations, bulk_stats.residual_norm);
end

% --- STAGE 2: NON-ISOBARIC CONFINED STEADY SWEEP ---
if isinf(r_cap)
    Pdew_Pa   = P_bulk_Pa;
    K_factors = K_bulk;
    Pcap_Pa   = 0.0;
    Pliq_Pa   = P_bulk_Pa;
    stats     = bulk_stats;
else
    fprintf('\nExecuting Stage 2: Tracing Confined Boundary at r_cap = %.2f nm [%s]...\n', r_nm, modelMode);
    
    if strcmpi(modelMode, 'FWI-Only')
        capModeSetting = 'PL';
    else
        capModeSetting = 'Pc';
    end
    
    % Thermodynamic Continuation Seeding
    if bulk_stats.converged && P_bulk_Pa > 0
        P_conf_seed = P_bulk_Pa;
        K_conf_seed = K_bulk;
    else
        P_conf_seed = P_guess_Pa;
        K_conf_seed = [];
    end
    
    [Pdew_Pa, K_factors, Pcap_Pa, Pliq_Pa, stats] = flash_conf.solveDewPoint(...
        T_K, P_conf_seed, z_feed, r_cap, ...
        'Solver', selectedSolver, ...
        'CapMode', capModeSetting, ...
        'PreconditionSS', useSSPrecond, ...
        'K_seed', K_conf_seed, ...
        'Pcap_seed', []);
end
execTime = toc;

%% 7. DIAGNOSTIC REPORT & THERMODYNAMIC SUMMARY
Pdew_MPa  = Pdew_Pa / 1e6;     Pdew_psi  = Pdew_Pa / psi2Pa;
Pbulk_MPa = P_bulk_Pa / 1e6;   Pbulk_psi = P_bulk_Pa / psi2Pa;
Pcap_MPa  = Pcap_Pa / 1e6;     Pcap_psi  = Pcap_Pa / psi2Pa;
Pliq_MPa  = Pliq_Pa / 1e6;     Pliq_psi  = Pliq_Pa / psi2Pa;
shift_psi = Pdew_psi - Pbulk_psi;
shift_MPa = Pdew_MPa - Pbulk_MPa;

fprintf('\n===================================================================\n');
fprintf('                 CONVERGED EQUILIBRIUM STATE REPORT                \n');
fprintf('===================================================================\n');
fprintf(' Physics Mode         : %s\n', modelMode);
fprintf(' Active Solver Engine : %s\n', selectedSolverName);
fprintf(' Solver Status        : %s (Completed in %.3f seconds)\n', ...
    upper(string(stats.converged)), execTime);
fprintf(' Total Iterations     : %d\n', stats.iterations);
fprintf(' Final Residual Norm  : %.2e\n', stats.residual_norm);
fprintf('-------------------------------------------------------------------\n');
fprintf(' Macroscopic Phase Boundaries:\n');
fprintf('   -> Bulk Dew Pressure (P_bulk) : %8.4f MPa  (%9.2f psia)\n', Pbulk_MPa, Pbulk_psi);
fprintf('   -> Confined Dew Point (P_dew) : %8.4f MPa  (%9.2f psia)\n', Pdew_MPa, Pdew_psi);
fprintf('   -> Liquid Phase Pressure (P_l): %8.4f MPa  (%9.2f psia)\n', Pliq_MPa, Pliq_psi);
fprintf('   -> Capillary Jump (P_cap)     : %8.4f MPa  (%9.2f psi)\n', Pcap_MPa, Pcap_psi);
fprintf('   -> Net Confined Shift (ΔPdew) : %+8.4f MPa  (%+9.2f psi)\n', shift_MPa, shift_psi);

if hasExpMatch
    fprintf('-------------------------------------------------------------------\n');
    fprintf(' Experimental Benchmark Comparison:\n');
    exp_bulk = expRow.("Bulk_Pdew");
    fprintf('   -> Experimental Bulk Pdew     : %9.2f psi\n', exp_bulk);
    if ismember('Pdew_Exp_lower', expRow.Properties.VariableNames)
        fprintf('   -> Experimental Confined Band : [%.1f, %.1f] psi\n', ...
            expRow.("Pdew_Exp_lower"), expRow.("Pdew_Exp_upper"));
        fprintf('   -> Experimental Shift Band    : [%+.1f, %+.1f] psi\n', ...
            expRow.("Pdew_Shift_Exp_lower"), expRow.("Pdew_Shift_Exp_upper"));
    end
end

fprintf('-------------------------------------------------------------------\n');
fprintf(' Component Equilibrium Split Factors (K_i = y_i / x_i):\n');
for i = 1:nc
    fprintf('   %-8s (z = %.4f)  ->  K = %12.6f  |  ln(K) = %9.4f\n', ...
        fluid.ComponentNames(i), z_feed(i), K_factors(i), log(max(K_factors(i), 1e-14)));
end
fprintf('===================================================================\n');
