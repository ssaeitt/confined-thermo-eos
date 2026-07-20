% =========================================================================
% CONFINED THERMODYNAMIC EOS MODELING: INTERACTIVE PDEW DRIVER
% Object-Oriented framework for VLE under nanoporous confinement
% =========================================================================
clc; clear; close all;
format shortG;

% 1. INITIALIZE WORKSPACE & NAMESPACE PATHS
if ~isfolder('src')
    error('ArchitectureError:MissingSource', 'The "src" directory was not found. Ensure you are running from the repository root.');
end
addpath('src');

% Locate database file across standard structural paths
xlsxPath = 'config/MixtureData.xlsx';
if ~isfile(xlsxPath)
    xlsxPath = 'MixtureData_2.xlsx'; % Fallback to local workspace root
    if ~isfile(xlsxPath)
        error('DataIOError:WorkbookNotFound', 'Thermodynamic database spreadsheet not found in config/ or root directory.');
    end
end

fprintf('===================================================================\n');
fprintf('   NANOPOROUS CONFINEMENT PDEW PREDICTION ENGINE (OOP V1.1)        \n');
fprintf('===================================================================\n\n');

%% 2. INTERACTIVE RESEARCH PARAMETER & MIXTURE CONFIGURATION

% --- Geological Formation Selection ---
fprintf('Select Target Geological Formation:\n');
fprintf('  1. Eagle Ford (EF2  - Carbonate-Rich, Low TOC)\n');
fprintf('  2. Barnett    (B1   - Silicate-Rich,  High TOC)\n');
rockChoice = input('Enter formation index (1-2): ');
rockMap = {'EF2', 'B1'};
if isempty(rockChoice) || ~ismember(rockChoice, [1, 2]), rockChoice = 1; end
selectedRock = rockMap{rockChoice};

% --- Hydrocarbon Mixture System Selection ---
fprintf('\nSelect Hydrocarbon Mixture System:\n');
fprintf('  1. C1-nC5  (Methane / n-Pentane)\n');
fprintf('  2. C1-nC8  (Methane / n-Octane)\n');
fprintf('  3. C1-nC10 (Methane / n-Decane)\n');
fprintf('  4. C1-nC5-nC10 (Ternary System)\n');
mixChoice = input('Enter mixture index (1-4): ');
mixMap = {'C1-nC5', 'C1-nC8', 'C1-nC10', 'C1-nC5-nC10'};
if isempty(mixChoice) || ~ismember(mixChoice, 1:4), mixChoice = 1; end
selectedMixName = mixMap{mixChoice};

% --- Thermodynamic State & Confinement Boundaries ---
T_C = input('\nEnter system temperature [°C] (default = 67.85): ');
if isempty(T_C), T_C = 67.85; end
T_K = T_C + 273.15;

r_nm = input('Enter pore radius [nm] (default = 10.0, enter Inf for bulk): ');
if isempty(r_nm), r_nm = 10.0; end
if isinf(r_nm)
    r_cap = Inf;
else
    r_cap = r_nm * 1e-9; % Convert nm to meters for SI thermodynamic consistency
end

P_guess_MPa = input('Enter initial Pdew guess [MPa] (default = 25.0): ');
if isempty(P_guess_MPa), P_guess_MPa = 25.0; end
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

compChoice = input(sprintf('Select feed composition index (1-%d): ', nRows + 1));
if isempty(compChoice) || compChoice < 1 || compChoice > (nRows + 1)
    compChoice = 1; % Default to first dataset row
end

if compChoice <= nRows
    rawStr = string(validRows.("mole_frac")(compChoice));
    % Strip brackets if necessary and evaluate numerical array safely
    cleanStr = strrep(strrep(rawStr, "[", ""), "]", "");
    z_feed = str2double(strsplit(cleanStr, ","));
    z_feed = z_feed(:) / sum(z_feed); % Enforce strict column orientation and normalization
    
    exp_Pdew_lower = validRows.("Pdew_Exp_lower")(compChoice);
    exp_Pdew_upper = validRows.("Pdew_Exp_upper")(compChoice);
    fprintf('      -> Selected Dataset Feed: z = [%s]\n', num2str(z_feed'));
else
    % Prompt user for custom vector
    fprintf('Enter molar composition as a vector matching [%s] (e.g., [0.85, 0.15]):\n', ...
        strjoin(fluid.ComponentNames, ', '));
    z_custom = input('z = ');
    if isempty(z_custom) || length(z_custom) ~= fluid.NC
        error('InputError:InvalidComposition', 'Custom vector must contain exactly %d elements.', fluid.NC);
    end
    z_feed = z_custom(:)' / sum(z_custom); % Enforce normalization
    exp_Pdew_lower = NaN; exp_Pdew_upper = NaN;
    fprintf('      -> Applied Custom Feed: z = [%s]\n', num2str(z_feed));
end
z_feed = z_feed(:); % Enforce column vector for matrix operations

%% 4. EXPLICIT TIP MATRIX (kijc) VIA INVERSE-DISPARITY CORRELATION
fprintf('[2/4] Computing Ternary Interaction Parameter (TIP) matrix...\n');

% 1. Define Correlation Fitting Coefficients (Inverse-Disparity Model)
A_corr = 90.5402886251168;
B_corr = 1.09588610509636;
C_corr = -150.825648642323;

nc = fluid.NC;
NA = 6.02214076e23; % Avogadro's Constant [1/mol]

% 2. Extract and Convert Energy Profiles to Molar Basis [J/mol]
eps_fluid_molar = fluid.LJ_Energy(:) * NA;
eps_rock_molar  = rock.EnergyVector(:);
w_norm          = rock.WeightVector(:);
w_TOC           = max(rock.w_TOC, 1e-4); % Guard against division-by-zero if TOC is 0

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
        val_corr = A_corr * E_star + (B_corr * E_star) / w_TOC + C_corr;
        
        kijc_matrix(r, c) = val_corr;
        kijc_matrix(c, r) = val_corr;
    end
end

fprintf('      -> TIP Matrix Configured (Inverse-Disparity Model | TOC = %.2f%%)\n', w_TOC * 100);

% 5. Instantiate Thermodynamic Engine with Computed Matrix
eos = thermo.ConfinedEOS(fluid, rock, 'kijc', kijc_matrix);

%% 5. INSTANTIATE NUMERICAL SOLVER ENGINES
fprintf('[3/4] Assembling non-isobaric Tangent Plane Distance (TPD) stability tester...\n');
stability = solvers.StabilityTester(eos, 'MaxIterations', 500, 'GradTolerance', 1e-6);

fprintf('[4/4] Initializing hybrid Newton-Raphson FlashEngine...\n');
flash = solvers.FlashEngine(eos, stability, 'MaxIterations', 1000, 'Tolerance', 1e-8);

%% 6. EXECUTE TWO-STAGE SATURATION POINT MOLECULAR CONTINUATION
tic;

% --- STAGE 1: UNCONFINED BULK CALIBRATION ANCHOR ---
fprintf('\nExecuting Stage 1: Solving Unconfined Bulk Saturation Boundary...\n');
[P_bulk, K_bulk, ~, ~, bulk_stats] = flash.solveDewPoint(...
    T_K, P_guess_Pa, z_feed, Inf, 'Solver', 'newton');

if bulk_stats.converged
    fprintf('      -> Bulk Stage Converged at %.4f MPa (%.2f psia)\n', P_bulk/1e6, P_bulk/6894.757);
    K_seed = K_bulk;
    P_start = P_bulk;
else
    fprintf('      -> Bulk Stage Stalled. Falling back to empirical seeds.\n');
    K_seed = [];
    P_start = P_guess_Pa;
end

% --- STAGE 2: NON-ISOBARIC CONFINED STEADY SWEEP ---
if isinf(r_cap)
    Pdew_Pa = P_bulk; K_factors = K_bulk; Pcap_Pa = 0.0; Pliq_Pa = P_bulk; stats = bulk_stats;
else
    fprintf('\nExecuting Stage 2: Tracing Confined Boundary at r_cap = %.2f nm...\n', r_nm);
    % Omit K_seed to force Stage 2 to use the internal StabilityTester
    % to find the true non-trivial confined phase coordinates at P_start
    [Pdew_Pa, K_factors, Pcap_Pa, Pliq_Pa, stats] = flash.solveDewPoint(...
        T_K, P_start, z_feed, r_cap, ...
        'Solver', 'newton', 'CapMode', 'Pc', ...
        'K_seed', [], 'Pcap_seed', 1e4);
end

execTime = toc;

%% 7. DIAGNOSTIC REPORT & THERMODYNAMIC SUMMARY
psi2Pa = 6894.757;
Pdew_MPa = Pdew_Pa / 1e6;   Pdew_psi = Pdew_Pa / psi2Pa;
Pcap_MPa = Pcap_Pa / 1e6;   Pcap_psi = Pcap_Pa / psi2Pa;
Pliq_MPa = Pliq_Pa / 1e6;   Pliq_psi = Pliq_Pa / psi2Pa;

fprintf('\n===================================================================\n');
fprintf('                 CONVERGED EQUILIBRIUM STATE REPORT                \n');
fprintf('===================================================================\n');
fprintf(' Solver Status      : %s (Completed in %.3f seconds)\n', upper(string(stats.converged)), execTime);
fprintf(' Total Iterations   : %d\n', stats.iterations);
fprintf(' Final Residual Norm: %.2e\n', stats.residual_norm);
fprintf('-------------------------------------------------------------------\n');
fprintf(' Macroscopic Phase Pressures:\n');
fprintf('   -> Vapor Dew Pressure (P_dew) : %8.4f MPa  (%9.2f psia)\n', Pdew_MPa, Pdew_psi);
fprintf('   -> Liquid Phase Pressure (P_l): %8.4f MPa  (%9.2f psia)\n', Pliq_MPa, Pliq_Pa / psi2Pa);
fprintf('   -> Capillary Discontinuity(Pc): %8.4f MPa  (%9.2f psi)\n', Pcap_MPa, Pcap_Pa / psi2Pa);
fprintf('-------------------------------------------------------------------\n');
fprintf(' Component Equilibrium Split Factors (K_i = y_i / x_i):\n');
for i = 1:nc
    fprintf('   %-8s (z = %.4f)  ->  K = %12.6f  |  ln(K) = %9.4f\n', ...
        fluid.ComponentNames(i), z_feed(i), K_factors(i), log(K_factors(i)));
end
fprintf('===================================================================\n');