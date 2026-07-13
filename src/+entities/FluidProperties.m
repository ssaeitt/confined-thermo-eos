classdef FluidProperties
    % FLUIDPROPERTIES Domain entity encapsulating pure component constants,
    % volume translation parameters, Lennard-Jones parameters, and the
    % fluid-fluid binary interaction matrix (BIP) for an NC-component mixture.
    
    properties (SetAccess = private)
        ComponentNames (1,:) string   % Component identifiers (e.g., ["C1", "nC5"])
        NC (1,1) double {mustBeInteger, mustBePositive} = 1 % Number of components with Fixed implicit default singularity guard
        
        % Core Cubic EOS Parameters
        Tc (1,:) double {mustBePositive, mustBeReal}    % Critical Temperature [K]
        Pc (1,:) double {mustBePositive, mustBeReal}    % Critical Pressure [Pa]
        omega (1,:) double {mustBeReal}                 % Acentric factor
        MW (1,:) double {mustBePositive, mustBeReal}    % Molecular Weight [g/mol]
        
        % Liquid Volume and Interfacial Tension Parameters
        Parachor (1,:) double {mustBePositive, mustBeReal}      % MacLeod-Sugden Parachor
        VolumeShift (1,:) double {mustBeReal}                   % Peneloux shift parameter (SE)
        
        % Microscopic Lennard-Jones Parameters for Fluid-Wall Alignment
        LJ_Size (1,:) double {mustBeNonnegative, mustBeReal}    % Collision diameter (sigma) [m]
        LJ_Energy (1,:) double {mustBeNonnegative, mustBeReal}  % Well depth potential (epsilon) [J]
        
        % Fluid-Fluid Cross Interactions
        BIPMatrix (:,:) double {mustBeReal}                     % Symmetric BIP matrix [NC x NC]
    end
    
    methods
        function obj = FluidProperties(names, tc, pc, omega, mw, parachor, vshift, lj_size, lj_energy, bip)
            % Constructor to explicitly map, size-check, and enforce array alignment
            if nargin > 0
                obj.NC = length(names);
                obj.ComponentNames = string(names);
                
                % Standardize all inputs to strict horizontal vectors
                obj.Tc = reshape(tc, 1, []);
                obj.Pc = reshape(pc, 1, []);
                obj.omega = reshape(omega, 1, []);
                obj.MW = reshape(mw, 1, []);
                obj.Parachor = reshape(parachor, 1, []);
                obj.VolumeShift = reshape(vshift, 1, []);
                obj.LJ_Size = reshape(lj_size, 1, []);
                obj.LJ_Energy = reshape(lj_energy, 1, []);
                
                % Assert and assign symmetric interaction bounds
                if size(bip, 1) ~= obj.NC || size(bip, 2) ~= obj.NC
                    error('FluidProperties:BIPDimensionMismatch', ...
                        'BIP Matrix size [%dx%d] must match the number of components NC = %d.', ...
                        size(bip, 1), size(bip, 2), obj.NC);
                end
                obj.BIPMatrix = bip;
                
                % Structural validation across all component indices
                obj.validateDimensions();
            end
        end
    end
    
    methods (Access = private)
        function validateDimensions(obj)
            % Enforces strict vector length compliance to prevent silent matrix execution crashes
            expectedLength = obj.NC;
            propertiesToCheck = {obj.Tc, obj.Pc, obj.omega, obj.MW, ...
                                 obj.Parachor, obj.VolumeShift, ...
                                 obj.LJ_Size, obj.LJ_Energy};
            for i = 1:length(propertiesToCheck)
                if length(propertiesToCheck{i}) ~= expectedLength
                    error('FluidProperties:VectorSizingFailure', ...
                        'Array dimensionality compromise detected. All characterization vectors must contain exactly NC = %d indices.', expectedLength);
                end
            end
        end
    end
    
    methods (Static)
        function obj = loadFromWorkbook(filePath, selectedMixName)
            % Automates ingestion and extraction across distinct worksheet layers
            
            % 1. Ingest Pure Component Properties
            optsProps = detectImportOptions(filePath, 'Sheet', 'Comp_Props');
            optsProps.VariableNamingRule = 'preserve';
            propsTable = readtable(filePath, optsProps);
            
            % If a specific mixture name is requested, filter rows to maintain explicit subset sequence
            if nargin > 1 && ~isempty(selectedMixName)
                components = strsplit(string(selectedMixName), '-');
                components = string(components);
                
                rowIndices = zeros(1, length(components));
                for i = 1:length(components)
                    idx = find(string(propsTable.("comp")) == components(i), 1);
                    if isempty(idx)
                        error('FluidProperties:ComponentNotFound', ...
                            'Component "%s" specified in mixture "%s" was not found in the Comp_Props worksheet.', ...
                            components(i), selectedMixName);
                    end
                    rowIndices(i) = idx;
                end
                propsTable = propsTable(rowIndices, :);
            end
            
            names = propsTable.("comp");
            tc = propsTable.("Tc");
            pc = propsTable.("Pc");
            omg = propsTable.("omega");
            mw = propsTable.("Mw");
            parachor = propsTable.("parachor");
            vshift = propsTable.("SE");
            lj_size = propsTable.("LJ_Size");
            lj_energy = propsTable.("LJ_Energy");
            
            % 2. Ingest Cross-Species BIP Intersect
            optsBIP = detectImportOptions(filePath, 'Sheet', 'BIP');
            optsBIP.VariableNamingRule = 'preserve';
            bipTable = readtable(filePath, optsBIP);
            
            % Extract numeric values, ignoring the leading text identifier column
            bipRaw = table2array(bipTable(:, 2:end));
            
            % 3. Map Subset/Sequence Integrity
            % Ensures that the sequence of components in 'BIP' aligns perfectly with the filtered 'names'
            bipComponentOrder = string(bipTable.Properties.VariableNames(2:end));
            
            sortIdx = zeros(1, length(names));
            for i = 1:length(names)
                idx = find(bipComponentOrder == string(names(i)), 1);
                if isempty(idx)
                    error('FluidProperties:BIPComponentNotFound', ...
                        'Component "%s" was not found in the BIP worksheet headers.', names(i));
                end
                sortIdx(i) = idx;
            end
            bipCorrected = bipRaw(sortIdx, sortIdx);
            
            % Instantiate the validated OOP Data Object
            obj = entities.FluidProperties(names, tc, pc, omg, mw, ...
                parachor, vshift, lj_size, lj_energy, bipCorrected);
        end
    end
end