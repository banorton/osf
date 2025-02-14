classdef ParaxialSystem
    properties
        numPlanes
        RIs
        focals
        marginalRay
        chiefRay
        distances
        imgPlaneDist
        mTransverse
        mAngular
        objHeight
        apStopDist
        apStopHeight
        apStopSurf
        apStopRatios
        CAs
    end
    
    methods
        function obj = ParaxialSystem(sim)
            % Construct a ParaxialSystem from a Sim object.
            %
            % Assumptions:
            %   - The refractive index is always 1.
            %   - For lens elements, use element.focalLength; for non-lenses, use inf.
            %   - sim.distances is a vector giving the distance in front of each element.
            %   - An object plane is added at the beginning and an image plane at the end.
            
            numElements = length(sim.elements);
            obj.numPlanes = numElements;
            
            % Refractive indices: assume air = 1 for all planes
            obj.RIs = ones(1, obj.numPlanes);
            
            % Extract focal lengths: for lenses use the given value, for others use inf.
            focals = zeros(1, numElements);
            for i = 1:numElements
                element = sim.elements{i};
                if strcmpi(element.elementType, 'lens')
                    focals(i) = element.focalLength;
                else
                    focals(i) = inf;
                end
            end
            obj.focals = focals;
            
            % Distances: add a zero at the beginning (object plane) and at the end (image plane)
            obj.distances = sim.distances;
            
            % Clear apertures: use a default value (e.g., 0.05) for all elements (object & image are infinite)
            obj.CAs = [repmat(0.05, 1, numElements)];
            
            % Object height: if defined in sim, use it; otherwise default to 0.02.
            if isprop(sim, 'objHeight')
                obj.objHeight = sim.objHeight;
            else
                obj.objHeight = 0.04;
            end
            
            % Other properties:
            obj.apStopDist = 0;
            obj.apStopHeight = 0;
            obj.apStopSurf = 0;
            
            % Initialize marginal and chief rays.
            obj.marginalRay.heights = zeros(1, obj.numPlanes);
            obj.marginalRay.angles  = zeros(1, obj.numPlanes);
            obj.marginalRay.heights(1) = 0;
            obj.marginalRay.angles(1)  = 0.05;
            
            obj.chiefRay.heights = zeros(1, obj.numPlanes);
            obj.chiefRay.angles  = zeros(1, obj.numPlanes);
            obj.chiefRay.heights(1) = obj.objHeight;
            obj.chiefRay.angles(1)  = 0;
        end
    end
    
    methods (Static)
        function system = solveSystem(system)
            system = ParaxialSystem.solveMarginal(system);
            system = ParaxialSystem.solveChief(system);
        end
        
        function system = solveMarginal(system)
            for i = 2:system.numPlanes
                system.marginalRay.heights(i) = system.marginalRay.heights(i-1) + system.marginalRay.angles(i-1) * system.distances(i);
                system.marginalRay.angles(i) = ((-system.marginalRay.heights(i) / system.focals(i)) + (system.RIs(i-1) * system.marginalRay.angles(i-1))) / system.RIs(i);
            end
            system.marginalRay.angles(end) = ((-system.marginalRay.heights(end) / system.focals(end)) + (system.RIs(end-1) * system.marginalRay.angles(end-1))) / system.RIs(end);
            system.imgPlaneDist = -system.marginalRay.heights(end-1) / system.marginalRay.angles(end-1);
            % system.distances(end) = system.imgPlaneDist;
            
            system.apStopRatios = abs(system.marginalRay.heights ./ system.CAs);
            system.marginalRay.heights = system.marginalRay.heights * (1 / max(system.apStopRatios));
        end
        
        function system = solveChief(system)
            for i = 2:system.numPlanes
                system.chiefRay.heights(i) = system.chiefRay.heights(i-1) + system.chiefRay.angles(i-1) * system.distances(i);
                system.chiefRay.angles(i) = ((-system.chiefRay.heights(i) / system.focals(i)) + (system.RIs(i-1) * system.chiefRay.angles(i-1))) / system.RIs(i);
            end
            system.mTransverse = system.chiefRay.heights(end) / system.chiefRay.heights(1);
            system.mAngular = system.chiefRay.angles(end) / system.chiefRay.angles(1);
        end
        
        function printSystem(system)
            fprintf('MARGINAL RAY\n');
            fprintf('  Heights\t  Angles\n');
            disp([system.marginalRay.heights', system.marginalRay.angles']);
            fprintf('\n');
            fprintf('CHIEF RAY\n');
            fprintf('  Heights\t  Angles\n');
            disp([system.chiefRay.heights', system.chiefRay.angles']);
            fprintf('\n');
            fprintf('IMAGE PLANE:     %f\n', system.imgPlaneDist);
            fprintf('IMAGE HEIGHT:    %f\n', system.chiefRay.heights(end));
            fprintf('TRANSVERSE MAG.: %f\n', system.mTransverse);
            fprintf('ANGULAR MAG.:    %f\n', system.mAngular);
        end
        
    end
end