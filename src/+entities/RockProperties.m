classdef RockProperties
    % ROCKPROPERTIES Domain entity handling mineral group classifications,
    % dynamic mass-weighted average surface energies, and normalization 
    % profiles for nanoporous confinement boundaries.
    
    properties (SetAccess = private)
        RockName (1,1) string          % Name identifier of target sample (e.g., "B1", "EF2")
        Theta (1,1) double {mustBeReal, mustBeNonnegative} % Macroscopic contact angle [degrees]

        % Unnormalized Raw TOC Mass Fraction (Preserved specifically for empirical TIP correlations)
        RawTOC (1,1) double {mustBeReal, mustBeNonnegative}

        % Normalized Mineral Group Mass Fractions (Sum = 1.0)
        w_Silicates (1,1) double {mustBeReal, mustBeNonnegative}
        w_Carbonates (1,1) double {mustBeReal, mustBeNonnegative}
        w_Clays (1,1) double {mustBeReal, mustBeNonnegative}
        w_Others (1,1) double {mustBeReal, mustBeNonnegative}
        w_TOC (1,1) double {mustBeReal, mustBeNonnegative}

        % Composite Group Surface Energy Parameters [J/mol or K-based basis]
        E_Silicates (1,1) double {mustBeReal, mustBeNonnegative}
        E_Carbonates (1,1) double {mustBeReal, mustBeNonnegative}
        E_Clays (1,1) double {mustBeReal, mustBeNonnegative}
        E_Others (1,1) double {mustBeReal, mustBeNonnegative}
        E_TOC (1,1) double {mustBeReal, mustBeNonnegative}
    end
    
    properties (Constant, Access = private)
        R double = 8.3144621            % Universal Gas Constant for K -> J/mol conversion
    end

    properties (Dependent)
        WeightVector (1,5) double      % Consolidated vector [Sil, Carb, Clay, Oth, TOC]
        EnergyVector (1,5) double      % Consolidated vector [Sil, Carb, Clay, Oth, TOC]
    end
    
    methods
        function obj = RockProperties(name, theta, weights, energies, rawTOC)
            % Constructor enforcing core property instantiation and validation
            if nargin > 0
                obj.RockName = string(name);
                obj.Theta = theta;

                if nargin > 4
                    obj.RawTOC = rawTOC;
                else
                    obj.RawTOC = weights(5); % Fallback default if raw fraction is omitted
                end
                
                % Unpack group configuration arrays
                obj.w_Silicates  = weights(1);
                obj.w_Carbonates = weights(2);
                obj.w_Clays      = weights(3);
                obj.w_Others     = weights(4);
                obj.w_TOC        = weights(5);
                
                obj.E_Silicates  = energies(1);
                obj.E_Carbonates = energies(2);
                obj.E_Clays      = energies(3);
                obj.E_Others     = energies(4);
                obj.E_TOC        = energies(5);
                
                % Assert mathematical mass normalization compliance
                totalMass = sum(weights);
                if abs(totalMass - 1.0) > 1e-6
                    error('RockProperties:NormalizationFailure', ...
                        'Constructed mineral mass sum (%f) must be normalized to exactly 1.0.', totalMass);
                end
            end
        end
        
        % Dependent Property Getters for Vectorized Matrix Operations Downstream
        function val = get.WeightVector(obj)
            val = [obj.w_Silicates, obj.w_Carbonates, obj.w_Clays, obj.w_Others, obj.w_TOC];
        end
        
        function val = get.EnergyVector(obj)
            val = [obj.E_Silicates, obj.E_Carbonates, obj.E_Clays, obj.E_Others, obj.E_TOC];
        end
    end
    
    methods (Static)
        function obj = loadFromWorkbook(filePath, targetRockName)
            % Ingests raw matrix allocations and applies grouping math rules dynamically
            opts = detectImportOptions(filePath, 'Sheet', 'Rock_Props');
            opts.VariableNamingRule = 'preserve';
            rawTable = readtable(filePath, opts);
            
            % Locate core structural definitions from row signatures
            rowTypes = string(rawTable.("Mineral"));
            energyIdx = find(rowTypes == "mineral_E", 1);
            weightRows = find(rowTypes == "mineral_w");
            
            if isempty(energyIdx) || isempty(weightRows)
                error('RockProperties:InvalidSheetFormat', 'Target sheet must contain explicit row identifiers.');
            end
            
            % Locate target row index corresponding to the specified rock sample formation
            targetRowIdx = [];
            for idx = reshape(weightRows, 1, [])
                if string(rawTable.("Rock")(idx)) == string(targetRockName)
                    targetRowIdx = idx;
                    break;
                end
            end
            
            if isempty(targetRowIdx)
                error('RockProperties:TargetRockNotFound', 'Rock sample matrix "%s" is absent from database entries.', targetRockName);
            end
            
            % Extract the macroscopic contact angle parameter boundary
            thetaVal = double(rawTable.("theta")(targetRowIdx));
            
            % Map Explicit Structural Columns to Structural Names
            minerals = ["Quartz", "Plagioclase", "K-feldspar", "Calcite", "Dolomite", ...
                        "Siderite", "Mica", "Kaolinite", "Chlorite", "Smectite", ...
                        "Pyrite", "Fluorapatite", "Gypsum", "TOC"];
                    
            E_raw = zeros(1, length(minerals));
            w_raw = zeros(1, length(minerals));
            
            for m = 1:length(minerals)
                E_raw(m) = double(rawTable.(minerals(m))(energyIdx));
                w_raw(m) = double(rawTable.(minerals(m))(targetRowIdx));
            end
            
            % Isolate Array Vectors Based on Structural Group Mappings
            silicate_indices  = [1, 2, 3];     % Quartz, Plagioclase, K-feldspar
            carbonate_indices = [4, 5, 6];     % Calcite, Dolomite, Siderite
            clay_indices      = [7, 8, 9, 10]; % Mica, Kaolinite, Chlorite, Smectite
            others_indices    = [11, 12, 13];  % Pyrite, Fluorapatite, Gypsum
            toc_index         = 14;            % Total Organic Carbon Explicit Group
            
            % Compute Raw Mass Accumulations per Allocation Block
            w_Sil_raw  = sum(w_raw(silicate_indices));
            w_Car_raw  = sum(w_raw(carbonate_indices));
            w_Clay_raw = sum(w_raw(clay_indices));
            w_Oth_raw  = sum(w_raw(others_indices));
            w_TOC_raw  = w_raw(toc_index);
            
            % Apply Group Energy Math Rule via Fully-Qualified Package Call
            E_Sil  = entities.RockProperties.calculateWeightedEnergy(w_raw(silicate_indices),  E_raw(silicate_indices));
            E_Car  = entities.RockProperties.calculateWeightedEnergy(w_raw(carbonate_indices), E_raw(carbonate_indices));
            E_Clay = entities.RockProperties.calculateWeightedEnergy(w_raw(clay_indices),      E_raw(clay_indices));
            E_Oth  = entities.RockProperties.calculateWeightedEnergy(w_raw(others_indices),    E_raw(others_indices));
            E_TOC  = E_raw(toc_index); 
            
            % Enforce Structural Normalization Across Collective Blocks
            sumRawWeights = w_Sil_raw + w_Car_raw + w_Clay_raw + w_Oth_raw + w_TOC_raw;
            if sumRawWeights <= 0
                error('RockProperties:ZeroMassError', 'Total evaluated rock matrix mass cannot be zero.');
            end
            
            normalizedWeights = [w_Sil_raw, w_Car_raw, w_Clay_raw, w_Oth_raw, w_TOC_raw] ./ sumRawWeights;

            % Apply unit conversion (Kelvin -> J/mol) by scaling with Gas Constant R
            R_val = entities.RockProperties.R;
            groupEnergies     = [E_Sil, E_Car, E_Clay, E_Oth, E_TOC] * R_val;
            
            % Return configured and normalized class instance
            obj = entities.RockProperties(targetRockName, thetaVal, normalizedWeights, groupEnergies, w_TOC_raw);
        end
        
        function E_avg = calculateWeightedEnergy(weights, energies)
            % Evaluates mass-weighted group energy profiles with division-by-zero protection
            sumW = sum(weights);
            if sumW > 0
                E_avg = sum(weights .* energies) / sumW;
            else
                E_avg = 0.0; % Handle edge cases where a formation lacks a specific group entirely
            end
        end
    end
end