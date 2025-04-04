classdef Sim < handle
    properties
        elements          % Array of Element objects
        distances         % Array of distances between elements
        resolution        % Size of each pixel (meters)
        fieldLength       % Size of the field of view (meters)
        samples           % Number of samples (always scalar)
        dim               % 1 or 2 - Dimensionality of the system
        wavelength            % Wavelength of light in meters (default 632 nm)
        paddingRatio      % Ratio of padding to add during propagation
    end

    methods

        %% CONSTRUCTOR
        function obj = Sim(resolution, fieldLength, varargin)
            p = inputParser;
            addRequired(p, 'resolution', @isnumeric);
            addRequired(p, 'fieldLength', @isnumeric);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'wavelength', 632e-9, @isnumeric);
            addParameter(p, 'paddingRatio', 1, @isnumeric);
            parse(p, resolution, fieldLength, varargin{:});

            obj.resolution = p.Results.resolution;
            obj.fieldLength = p.Results.fieldLength;
            obj.dim = p.Results.dim;
            obj.wavelength = p.Results.wavelength;
            obj.paddingRatio = p.Results.paddingRatio;

            obj.samples = round(obj.fieldLength / obj.resolution);
            obj.elements = {};
            obj.distances = [];
        end

        %% ELEMENT MANAGEMENT

        function addElement(obj, dist, element)
            % Add an optical element to the system with a specific distance
            % from the previous element.

            if nargin < 3
                dist = 0;
            end

            if obj.dim ~= element.dim
                error('Dimensionality of element must match the system dimensionality.');
            end

            obj.elements{end+1} = element;
            obj.distances(end+1) = dist;
        end

        function lens = addLens(obj, dist, focalLength, varargin)
            lens = osf.elements.Lens(focalLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, lens);
        end

        function diffuser = addDiffuser(obj, dist, roughness, correlationLength, varargin)
            diffuser = osf.elements.Diffuser(roughness, correlationLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, diffuser);
        end

        function aperture = addAperture(obj, dist, varargin)
            aperture = osf.elements.Aperture('dim', obj.dim, varargin{:});
            obj.addElement(dist, aperture);
        end

        function plane = addPlane(obj, dist, varargin)
            plane = osf.elements.Plane('dim', obj.dim, varargin{:});
            obj.addElement(dist, plane);
        end

        function filter = addFilter(obj, dist, field, varargin)
            filter = osf.elements.Filter(field, varargin{:});
            obj.addElement(dist, filter);
        end

        function grating = addGrating(obj, dist, linesPerMM, varargin)
            grating = osf.elements.Grating(linesPerMM, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, grating);
        end

        function source = addSource(obj, varargin)
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
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'resolution', [obj.samples, obj.samples]);
            addParameter(p, 'pixelPitch', obj.resolution);
            parse(p, varargin{:});

            detector = osf.elements.Detector(p.Results.resolution, p.Results.pixelPitch, 'dim', obj.dim, p.Unmatched);
            obj.addElement(dist, detector);
        end

        %% PROPAGATION METHODS
        % These methods are typically wrappers around the propagation calculation methods.

        function result = propToIndex(obj, field, targetIndex, varargin)
            p = inputParser;
            addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            addParameter(p, 'collect', false, @islogical);
            parse(p, varargin{:});

            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;
            collect = p.Results.collect;

            if targetIndex > length(obj.elements)
                error('Target index exceeds the number of elements in the system.');
            end

            if verbose
                fprintf('Starting propagation through all elements in the system\n');
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

                    fprintf('Propagated to element %d (%s) at distance: %.3e m\n', i, elementTitle, cumulativeDist);
                    field.show(title=sprintf('%s\nDist: %.0f mm', elementTitle, 1000 * cumulativeDist), figPosition=[750 200 450 700]);
                end
            end

            if collect
                result = collectedFields;
            else
                result = field;
            end
        end

        function result = prop(obj, field, varargin)
            lastElementIndex = length(obj.elements);
            if lastElementIndex > 0
                result = obj.propToIndex(field, lastElementIndex, varargin{:});
            end
        end

        function [field, currDist] = propToDist(obj, field, targetDist, varargin)
            % Propagate to a specific distance in the system

            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'Field'));
            addRequired(p, 'targetDist', @isnumeric);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, targetDist, varargin{:});

            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;

            currDist = 0;

            if verbose
                fprintf('Starting propagation from distance 0 to target distance: %.3e m\n', targetDist);
            end

            if isempty(obj.elements)
                % If no elements, propagate directly to target distance
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
                        fprintf('Propagated to distance: %.3e m (no elements present)\n', currDist);
                        % field.disp('title', 'Field after propagation (no elements)');
                        field.show('title', 'Field after propagation (no elements)');
                    end
                end
                return;
            end

            % Determine the last element index before the target distance
            cumulativeDist = cumsum(obj.distances);
            lastElementIndex = find(cumulativeDist <= targetDist, 1, 'last');

            % Propagate to last element before target distance
            if ~isempty(lastElementIndex)
                field = obj.propToIndex(field, lastElementIndex, 'verbose', verbose, 'propMethod', propMethod);
                currDist = cumulativeDist(lastElementIndex);
            end

            % Propagate the remaining distance to reach the target
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
                    fprintf('Propagated remaining distance to target: %.3e m\n', currDist);
                    % field.disp('title', 'Field after propagation to target distance');
                    field.show('title', 'Field after propagation to target distance');
                end
            end
        end

        function field = propToElement(obj, field, targetName, varargin)
            % Propagate to an element by its name

            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'osf.Field'));
            addRequired(p, 'targetName', @ischar);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, targetName, varargin{:});
            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;

            elementNames = cellfun(@(e) e.name, obj.elements, 'UniformOutput', false);
            elementIndex = find(strcmp(elementNames, targetName), 1);

            if isempty(elementIndex)
                error('Element with the specified name not found.');
            end

            if verbose
                fprintf('Starting propagation to element named: %s\n', targetName);
            end

            field = obj.propToIndex(field, elementIndex, 'verbose', verbose, 'propMethod', propMethod);

            if verbose
                cumulativeDist = sum(obj.distances(1:elementIndex));
                fprintf('Finished propagating to element named: %s at distance: %.3e m\n', targetName, cumulativeDist);
                % field.disp('title', sprintf('Field after propagation to element (%s), Distance: %.3e m', targetName, cumulativeDist));
                field.show('title', sprintf('Field after propagation to element (%s), Distance: %.3e m', targetName, cumulativeDist));
            end
        end

        function img = propScan(obj, field, numSteps, varargin)
            p = inputParser;
            addParameter(p, 'live', false, @islogical);
            parse(p, varargin{:});
            live = p.Results.live;

            totalDistance = sum(obj.distances);
            dz = totalDistance / numSteps;
            % fprintf('Total propagation steps: %d (dz = %.3e m)\n', numSteps, dz);

            img = zeros(obj.samples, numSteps);
            currentDistance = 0;
            nextElementIdx = 1;

            if live
                f = figure('Position', [313 208 1315 603]);
            end

            for step = 1:numSteps
                % fprintf('Propagation step %d of %d\n', step, numSteps);

                intensity = abs(field.getComplexField()).^2;
                if obj.dim == 1
                    vec = intensity(:);
                else
                    vec = mean(intensity, 2);
                end

                vec = vec / max(vec(:) + eps);
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

        %% PROPAGATION CALCULATION METHODS

        function field = angularSpectrumPropagation(obj, field, z)
            if z == 0
                return;
            end

            % Get complex input field and add padding.
            u_out = obj.addPadding(field.getComplexField());

            % Parameters for propagation.
            [Nx, Ny] = size(u_out);
            dx = obj.resolution; dy = obj.resolution;
            wavelength = obj.wavelength;

            % Frequency axis computation.
            fx = (-Nx/2:Nx/2-1) / (Nx * dx);
            fy = (-Ny/2:Ny/2-1) / (Ny * dy);
            [FX_sqr, FY_sqr] = meshgrid(fx.^2, fy.^2);

            % Compute longitudinal wavevector component k_z.
            k_z = 2 * pi * sqrt((1 / wavelength^2) - FX_sqr - FY_sqr);
            clear FX_sqr FY_sqr;

            % Compute transfer function.
            H = exp(1i * k_z * z);
            clear k_z;

            % Find the angular spectrum.
            u_out = fftshift(fft2(u_out));

            % Propagate the angular spectrum using the transfer function.
            u_out = u_out .* H;

            % Inverse Fourier transform to get the propagated field.
            u_out = obj.removePadding(ifft2(ifftshift(u_out)));

            field = field.setComplexField(u_out);
        end

        function field = angularSpectrumPropagation1D(obj, field, z)
            % Angular Spectrum Propagation (1D)
            dx = obj.resolution;
            uin = field.getComplexField();

            % Pad the input array
            uin_padded = obj.addPadding(uin);

            Nx = length(uin_padded);
            k = 2 * pi / obj.wavelength;

            dfx = 1 / (Nx * dx);
            fx = (-Nx / 2 : Nx / 2 - 1) * dfx;

            % Kernel for propagation in 1D
            kernel = exp(1i * k * z * sqrt(1 - (obj.wavelength * fx).^2));
            ftu = fftshift(fft(uin_padded));
            ftu = ftu .* kernel;
            uout_padded = ifft(ifftshift(ftu));

            uout = obj.removePadding(uout_padded);

            field = field.setComplexField(uout);
        end

        %% PADDING METHODS

        function paddedField = addPadding(obj, field)
            % Adds zero-padding to the input field based on paddingRatio
            paddingSize = round(obj.paddingRatio * obj.samples);
            if obj.dim == 1
                paddedField = padarray(field, [0, paddingSize], 0, 'both');
            else
                paddedField = padarray(field, [paddingSize, paddingSize], 0, 'both');
            end
        end

        function croppedField = removePadding(obj, paddedField)
            % Removes zero-padding from the input field
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

        %% UTILITY METHODS

        function obj = show(obj, varargin)
            osf.Dashboard([1,1], varargin{:}).show(obj);
        end

        function field = newField(obj, varargin)
            % Create a new Field instance with system dimensions
            % Optional input:
            %   'cmap' - colormap for the phase plot (inherits from Field if not specified)

            p = inputParser;
            addOptional(p, 'cmap', [], @(x) ischar(x) && ismember(x, ...
            {'turbo', 'parula', 'jet', 'hot', 'gray', 'cool', ...
            'spring', 'summer', 'autumn', 'winter', 'bone', 'copper', 'pink'}));
            parse(p, varargin{:});

            field = osf.Field(obj.dim, obj.fieldLength, obj.resolution, obj.wavelength);

            % Only override cmap if the user provides an input
            if ~isempty(p.Results.cmap)
                field.cmap = p.Results.cmap;
            end
        end

        function print(obj)
            %   Prints simulation parameters in a nicely formatted manner.
            %   Length measurements are converted to mm (or um for resolution)
            %   with no digits after the decimal, and element names are printed in a
            %   bracketed list. If an element's name is empty, its type is used instead.

            distances_mm   = obj.distances * 1e3;
            fieldLength_mm = obj.fieldLength * 1e3;
            resolution_um  = obj.resolution * 1e6;
            wavelength_nm      = obj.wavelength * 1e9;

            % Build a cell array of element names (or types if name is empty)
            numEl = numel(obj.elements);
            names = cell(1, numEl);
            for k = 1:numEl
                if isempty(obj.elements{k}.name)
                    names{k} = obj.elements{k}.elementType;
                else
                    names{k} = obj.elements{k}.name;
                end
            end
            % Join names with a comma and enclose in square brackets
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
