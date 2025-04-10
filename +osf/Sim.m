classdef Sim < handle
    % Sim class for managing optical simulations.
    % This class handles optical elements, field propagation, and field creation
    % for Fourier optics simulations.

    properties
        elements          % Cell array of optical element objects
        distances         % Array of distances (meters) between consecutive elements
        resolution        % Pixel size (meters)
        fieldLength       % Field of view (meters)
        samples           % Number of samples (scalar)
        dim               % Dimensionality (1 for 1D or 2 for 2D)
        wavelength        % Wavelength (meters; default 632e-9)
        paddingRatio      % Ratio of zero-padding added during propagation
    end

    methods

        %% Constructor
        function obj = Sim(resolution, fieldLength, varargin)
            % Constructs a new simulation.
            %
            % Required
            %   resolution : Pixel size in meters.
            %   fieldLength : Field of view in meters.
            %
            % Parameters (Name-Value Pairs)
            %   'dim' : Dimensionality (1 or 2; default: 2).
            %   'wavelength' : Wavelength in meters (default: 632e-9).
            %   'paddingRatio' : Padding ratio for propagation (default: 1).
            %
            % Example
            %   sim = Sim(5e-6, 1e-3, 'dim', 2);
            p = inputParser;
            addRequired(p, 'resolution', @isnumeric);
            addRequired(p, 'fieldLength', @isnumeric);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'wavelength', 632e-9, @isnumeric);
            addParameter(p, 'paddingRatio', 1, @isnumeric);
            parse(p, resolution, fieldLength, varargin{:});

            obj.resolution    = p.Results.resolution;
            obj.fieldLength   = p.Results.fieldLength;
            obj.dim           = p.Results.dim;
            obj.wavelength    = p.Results.wavelength;
            obj.paddingRatio  = p.Results.paddingRatio;

            obj.samples       = round(obj.fieldLength / obj.resolution);
            obj.elements      = {};
            obj.distances     = [];
        end

        %% Element Management

        function addElement(obj, dist, element)
            % Adds an optical element to the simulation.
            %
            % Required
            %   element : Optical element object to add.
            %
            % Optionals
            %   dist : Distance in meters from the previous element (default: 0).
            if nargin < 3
                dist = 0;
            end

            if obj.dim ~= element.dim
                error('Element dimensionality must match the simulation system.');
            end

            obj.elements{end+1} = element;
            obj.distances(end+1) = dist;
        end

        function removeElement(obj, identifier)
            % Removes an element by index or name from the simulation and updates lens numbering.
            %
            % Required
            %   identifier : Either an integer index of the element to remove or a string representing the element's name.

            if isnumeric(identifier)
                index = identifier;
            elseif ischar(identifier) || isstring(identifier)
                index = [];
                for i = 1:numel(obj.elements)
                    if isprop(obj.elements{i}, 'name') && strcmp(obj.elements{i}.name, char(identifier))
                        index = i;
                        break;
                    end
                end
                if isempty(index)
                    error('Element with name "%s" not found.', char(identifier));
                end
            else
                error('Identifier must be numeric (index) or a string (element name).');
            end

            if ~(isscalar(index) && isnumeric(index) && index == floor(index))
                error('Index must be an integer scalar.');
            end
            if index < 1 || index > numel(obj.elements)
                error('Invalid element index. Must be between 1 and %d.', numel(obj.elements));
            end

            obj.elements(index) = [];
            obj.distances(index) = [];

            lensCount = 0;
            for i = 1:numel(obj.elements)
                if isprop(obj.elements{i}, 'type') && strcmp(obj.elements{i}.type, 'lens')
                    lensCount = lensCount + 1;
                    obj.elements{i}.name = sprintf('Lens %d', lensCount);
                end
            end
        end

        function lens = addLens(obj, dist, focalLength, varargin)
            % Creates and adds a lens element.
            %
            % Required
            %   focalLength : Focal length of the lens in meters.
            %
            % Optionals
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for the lens.

            count = 0;
            for i = 1:numel(obj.elements)
                if isprop(obj.elements{i}, 'type') && strcmp(obj.elements{i}.type, 'lens')
                    count = count + 1;
                end
            end

            lensName = sprintf('Lens %d', count + 1);
            lens = osf.elements.Lens(focalLength, 'dim', obj.dim, 'name', lensName, varargin{:});
            obj.addElement(dist, lens);
        end

        function diffuser = addDiffuser(obj, dist, roughness, correlationLength, varargin)
            % Creates and adds a diffuser element.
            %
            % Required
            %   roughness : Surface roughness in meters.
            %   correlationLength : Correlation length in meters.
            %
            % Optionals
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for the diffuser.
            diffuser = osf.elements.Diffuser(roughness, correlationLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, diffuser);
        end

        function aperture = addAperture(obj, dist, varargin)
            % Creates and adds an aperture element.
            %
            % Required
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for the aperture.
            aperture = osf.elements.Aperture('dim', obj.dim, varargin{:});
            obj.addElement(dist, aperture);
        end

        function plane = addPlane(obj, dist, varargin)
            % Creates and adds a plane element (e.g., detector plane).
            %
            % Required
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs (e.g., 'name' for labeling).
            plane = osf.elements.Plane('dim', obj.dim, varargin{:});
            obj.addElement(dist, plane);
        end

        function filter = addFilter(obj, dist, field, varargin)
            % Creates and adds a filter element.
            %
            % Required
            %   dist : Distance in meters from the previous element.
            %   field : Field object defining the filter characteristics.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for filter configuration.
            filter = osf.elements.Filter(field, varargin{:});
            obj.addElement(dist, filter);
        end

        function grating = addGrating(obj, dist, linesPerMM, varargin)
            % Creates and adds a grating element.
            %
            % Required
            %   linesPerMM : Grating spatial frequency in lines per mm.
            %
            % Optionals
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for the grating.
            grating = osf.elements.Grating(linesPerMM, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, grating);
        end

        function source = addSource(obj, varargin)
            % Creates and adds a source element.
            %
            % Optionals
            %   wavelength : Wavelength in meters. If not provided, defaults to obj.wavelength.
            %
            % Parameters (Name-Value Pairs)
            %   Additional name-value pairs for the source.
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'wavelength', []);
            parse(p, varargin{:});

            if isempty(p.Results.wavelength)
                wavelength = obj.wavelength;
            else
                wavelength = p.Results.wavelength;
                obj.wavelength = wavelength;
            end

            source = osf.elements.Source(wavelength, 'dim', obj.dim, p.Unmatched);
            obj.addElement(0, source);
        end

        function detector = addDetector(obj, dist, varargin)
            % Creates and adds a detector element.
            %
            % Required
            %   dist : Distance in meters from the previous element.
            %
            % Parameters (Name-Value Pairs)
            %   'resolution' : Resolution as [Nx, Ny] (default: [samples, samples]).
            %   'pixelPitch' : Size of a pixel (default: simulation resolution).
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'resolution', [obj.samples, obj.samples]);
            addParameter(p, 'pixelPitch', obj.resolution);
            parse(p, varargin{:});

            detector = osf.elements.Detector(p.Results.resolution, p.Results.pixelPitch, 'dim', obj.dim, p.Unmatched);
            obj.addElement(dist, detector);
        end

        %% Propagation Methods

        function result = propToIndex(obj, field, targetIndex, varargin)
            % Propagates a field to a specified element index.
            %
            % Required
            %   field : Field object to propagate.
            %   targetIndex : Target element index (integer).
            %
            % Optionals
            %   'verbose' : Logical flag for verbose output (default: false).
            %   'propMethod' : Propagation method ('as' for angular spectrum or 'rs' for Rayleigh-Sommerfeld; default: 'as').
            %   'collect' : Logical flag to collect and return fields at each step (default: false).
            p = inputParser;
            addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            addParameter(p, 'collect', false, @islogical);
            parse(p, varargin{:});

            verbose    = p.Results.verbose;
            propMethod = p.Results.propMethod;
            collect    = p.Results.collect;

            if targetIndex > length(obj.elements)
                error('Target index exceeds number of elements.');
            end

            if verbose
                fprintf('Starting propagation through %d elements\n', targetIndex);
            end

            cumulativeDist = 0;
            collectedFields = {};

            for i = 1:targetIndex
                segmentDist = obj.distances(i);
                if obj.dim == 1
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation1D(field, segmentDist);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation1D(field, segmentDist);
                    end
                else
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation(field, segmentDist);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation(field, segmentDist);
                    end
                end
                cumulativeDist = cumulativeDist + segmentDist;
                field = obj.elements{i}.apply(field);

                if collect
                    collectedFields{end+1} = field;
                end

                if verbose
                    elementName = obj.elements{i}.name;
                    if isempty(elementName)
                        elementTitle = sprintf('Element %d', i);
                    else
                        elementTitle = elementName;
                    end
                    fprintf('Propagated to %s at %.3e m\n', elementTitle, cumulativeDist);
                    field.show(title=sprintf('%s\nDistance: %.0f mm', elementTitle, 1000 * cumulativeDist), ...
                    figPosition=[750 200 450 700]);
                end
            end

            if collect
                result = collectedFields;
            else
                result = field;
            end
        end

        function result = prop(obj, field, varargin)
            % Propagates a field through all elements.
            %
            % Required
            %   field : Field object to propagate.
            %
            % Optionals
            %   All optional parameters are forwarded to propToIndex.
            lastElementIndex = length(obj.elements);
            if lastElementIndex > 0
                result = obj.propToIndex(field, lastElementIndex, varargin{:});
            end
        end

        function [field, currDist] = propToDist(obj, field, targetDist, varargin)
            % Propagates a field to a specific distance.
            %
            % Required
            %   field : Field object to propagate.
            %   targetDist : Target propagation distance in meters.
            %
            % Optionals
            %   'verbose' : Logical flag for verbose output (default: false).
            %   'propMethod' : Propagation method ('as' or 'rs'; default: 'as').
            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'Field'));
            addRequired(p, 'targetDist', @isnumeric);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, targetDist, varargin{:});

            verbose    = p.Results.verbose;
            propMethod = p.Results.propMethod;
            currDist   = 0;

            if verbose
                fprintf('Propagating from 0 to %.3e m\n', targetDist);
            end

            if isempty(obj.elements)
                if targetDist > 0
                    if obj.dim == 1
                        if strcmp(propMethod, 'as')
                            field = obj.angularSpectrumPropagation1D(field, targetDist);
                        elseif strcmp(propMethod, 'rs')
                            field = obj.rayleighSommerfeldPropagation1D(field, targetDist);
                        end
                    else
                        if strcmp(propMethod, 'as')
                            field = obj.angularSpectrumPropagation(field, targetDist);
                        elseif strcmp(propMethod, 'rs')
                            field = obj.rayleighSommerfeldPropagation(field, targetDist);
                        end
                    end
                    currDist = targetDist;
                    if verbose
                        fprintf('Propagated %.3e m (no elements present)\n', currDist);
                        field.show('title', 'Field after propagation (no elements)');
                    end
                end
                return;
            end

            cumulativeDist = cumsum(obj.distances);
            lastElementIndex = find(cumulativeDist <= targetDist, 1, 'last');

            if ~isempty(lastElementIndex)
                field = obj.propToIndex(field, lastElementIndex, 'verbose', verbose, 'propMethod', propMethod);
                currDist = cumulativeDist(lastElementIndex);
            end

            if targetDist > currDist
                remainingDist = targetDist - currDist;
                if obj.dim == 1
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation1D(field, remainingDist);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation1D(field, remainingDist);
                    end
                else
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation(field, remainingDist);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation(field, remainingDist);
                    end
                end
                currDist = targetDist;
                if verbose
                    fprintf('Propagated remaining %.3e m to reach target\n', remainingDist);
                    field.show('title', 'Field after reaching target distance');
                end
            end
        end

        function field = propToElement(obj, field, targetName, varargin)
            % Propagates a field to a specified element by name.
            %
            % Required
            %   field : Field object to propagate.
            %   targetName : Name of the target element (string).
            %
            % Optionals
            %   'verbose' : Logical flag for verbose output (default: false).
            %   'propMethod' : Propagation method ('as' or 'rs'; default: 'as').
            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'osf.Field'));
            addRequired(p, 'targetName', @ischar);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, targetName, varargin{:});
            verbose    = p.Results.verbose;
            propMethod = p.Results.propMethod;

            elementNames = cellfun(@(e) e.name, obj.elements, 'UniformOutput', false);
            elementIndex = find(strcmp(elementNames, targetName), 1);
            if isempty(elementIndex)
                error('Element with name "%s" not found.', targetName);
            end

            if verbose
                fprintf('Propagating to element "%s"\n', targetName);
            end

            field = obj.propToIndex(field, elementIndex, 'verbose', verbose, 'propMethod', propMethod);

            if verbose
                cumulativeDist = sum(obj.distances(1:elementIndex));
                fprintf('Reached "%s" at %.3e m\n', targetName, cumulativeDist);
                field.show('title', sprintf('Field at "%s", Dist: %.3e m', targetName, cumulativeDist));
            end
        end

        function img = propScan(obj, field, numSteps, varargin)
            % Propagates a field stepwise and records its intensity.
            %
            % Required
            %   field : Field object to propagate.
            %   numSteps : Number of propagation steps.
            %
            % Optionals
            %   'live' : Logical flag for live updating during propagation (default: false).
            p = inputParser;
            addParameter(p, 'live', false, @islogical);
            parse(p, varargin{:});
            live = p.Results.live;

            totalDistance = sum(obj.distances);
            dz = totalDistance / numSteps;
            img = zeros(obj.samples, numSteps);
            currentDistance = 0;
            nextElementIdx = 1;

            if live
                f = figure('Position', [313 208 1315 603]);
            end

            for step = 1:numSteps
                intensity = abs(field.getComplexField()).^2;
                if obj.dim == 1
                    vec = intensity(:);
                else
                    vec = mean(intensity, 2);
                end
                vec = vec / (max(vec(:)) + eps);
                img(:, step) = vec;

                if obj.dim == 1
                    field = obj.angularSpectrumPropagation1D(field, dz);
                else
                    field = obj.angularSpectrumPropagation(field, dz);
                end
                currentDistance = currentDistance + dz;

                while nextElementIdx <= numel(obj.distances) && ...
                    currentDistance >= sum(obj.distances(1:nextElementIdx))
                    field = obj.elements{nextElementIdx}.apply(field);
                    nextElementIdx = nextElementIdx + 1;
                end

                if live
                    imagesc(img);
                    colormap turbo;
                    osf.vis.applyTheme(f);
                    axis off;
                    drawnow;
                end
            end
        end

        %% Propagation Calculation Methods

        function field = angularSpectrumPropagation(obj, field, z)
            % Propagates a field using the 2D angular spectrum method.
            %
            % Required
            %   field : Field object to propagate.
            %   z : Propagation distance in meters.
            if z == 0
                return;
            end

            u_out = obj.addPadding(field.getComplexField());
            [Nx, Ny] = size(u_out);
            dx = obj.resolution;
            dy = obj.resolution;
            wavelength = obj.wavelength;

            fx = (-Nx/2:Nx/2-1) / (Nx * dx);
            fy = (-Ny/2:Ny/2-1) / (Ny * dy);
            [FX_sqr, FY_sqr] = meshgrid(fx.^2, fy.^2);

            k_z = 2 * pi * sqrt((1 / wavelength^2) - FX_sqr - FY_sqr);
            clear FX_sqr FY_sqr;

            H = exp(1i * k_z * z);
            clear k_z;

            u_out = fftshift(fft2(u_out));
            u_out = u_out .* H;
            u_out = obj.removePadding(ifft2(ifftshift(u_out)));

            field = field.setComplexField(u_out);
        end

        function field = angularSpectrumPropagation1D(obj, field, z)
            % Propagates a field using the 1D angular spectrum method.
            %
            % Required
            %   field : Field object to propagate.
            %   z : Propagation distance in meters.
            dx = obj.resolution;
            uin = field.getComplexField();
            uin_padded = obj.addPadding(uin);
            Nx = length(uin_padded);
            k = 2 * pi / obj.wavelength;
            dfx = 1 / (Nx * dx);
            fx = (-Nx / 2 : Nx / 2 - 1) * dfx;
            kernel = exp(1i * k * z * sqrt(1 - (obj.wavelength * fx).^2));
            ftu = fftshift(fft(uin_padded));
            ftu = ftu .* kernel;
            uout_padded = ifft(ifftshift(ftu));
            uout = obj.removePadding(uout_padded);
            field = field.setComplexField(uout);
        end

        %% Padding Methods

        function paddedField = addPadding(obj, field)
            % Adds zero-padding to a field.
            %
            % Required
            %   field : Complex field matrix.
            %
            % Output
            %   paddedField : Field with zero-padding added.
            paddingSize = round(obj.paddingRatio * obj.samples);
            if obj.dim == 1
                paddedField = padarray(field, [0, paddingSize], 0, 'both');
            else
                paddedField = padarray(field, [paddingSize, paddingSize], 0, 'both');
            end
        end

        function croppedField = removePadding(obj, paddedField)
            % Removes zero-padding from a field.
            %
            % Required
            %   paddedField : Field matrix with padding.
            %
            % Output
            %   croppedField : Field matrix with padding removed.
            paddingSize = round(obj.paddingRatio * obj.samples);
            if obj.dim == 1
                startIndexX = paddingSize + 1;
                endIndexX = startIndexX + obj.samples - 1;
                croppedField = paddedField(startIndexX:endIndexX);
            else
                startIndexY = paddingSize + 1;
                endIndexY = startIndexY + obj.samples - 1;
                startIndexX = paddingSize + 1;
                endIndexX = startIndexX + obj.samples - 1;
                croppedField = paddedField(startIndexY:endIndexY, startIndexX:endIndexX);
            end
        end

        %% Utility Methods

        function obj = show(obj, varargin)
            % Displays simulation parameters via a dashboard.
            %
            % Optionals
            %   Additional name-value pairs for the dashboard display.
            osf.Dashboard([1,1], varargin{:}).show(obj);
        end

        function field = newField(obj, varargin)
            % Creates a new Field object.
            %
            % Optionals
            %   'cmap' : Colormap for phase display (e.g., 'turbo', 'parula', etc.).
            p = inputParser;
            addOptional(p, 'cmap', [], @(x) ischar(x) && ismember(x, ...
            {'turbo', 'parula', 'jet', 'hot', 'gray', 'cool', ...
            'spring', 'summer', 'autumn', 'winter', 'bone', 'copper', 'pink'}));
            parse(p, varargin{:});

            field = osf.Field(obj.dim, obj.fieldLength, obj.resolution, obj.wavelength);

            if ~isempty(p.Results.cmap)
                field.cmap = p.Results.cmap;
            end
        end

        function print(obj)
            % Displays a summary of the simulation parameters.
            %
            % Required
            %   No inputs.
            %
            % Output
            %   Prints simulation parameters including element list, distances,
            %   resolution, field length, number of samples, dimensionality,
            %   wavelength, and padding ratio.
            distances_mm   = obj.distances * 1e3;
            fieldLength_mm = obj.fieldLength * 1e3;
            resolution_um  = obj.resolution * 1e6;
            wavelength_nm  = obj.wavelength * 1e9;

            numEl = numel(obj.elements);
            names = cell(1, numEl);
            for k = 1:numEl
                if isempty(obj.elements{k}.name)
                    names{k} = obj.elements{k}.elementType;
                else
                    names{k} = obj.elements{k}.name;
                end
            end
            elementsStr = ['[', strjoin(names, ', '), ']'];

            fprintf('\nSim Parameters:\n');
            fprintf('-------------------------------\n');
            fprintf('  Elements:       %s\n', elementsStr);
            fprintf('  Distances:      %s\n', mat2str(distances_mm));
            fprintf('  Resolution:     %.2f um\n', resolution_um);
            fprintf('  Field Length:   %.1f mm\n', fieldLength_mm);
            fprintf('  Samples:        %d\n', obj.samples);
            fprintf('  Dimensionality: %d\n', obj.dim);
            fprintf('  Wavelength:     %d nm\n', wavelength_nm);
            fprintf('  Padding Ratio:  %.1f\n\n', obj.paddingRatio);
        end

    end
end
