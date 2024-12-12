classdef Sim < handle
    properties
        elements          % Array of Element objects
        distances         % Array of distances between elements
        resolution        % Size of each pixel (meters)
        fieldLength       % Size of the field of view (meters)
        samples           % Number of samples (always scalar)
        dim               % 1 or 2 - Dimensionality of the system
        lambda            % Wavelength of light in meters (default 632 nm)
        paddingRatio      % Ratio of padding to add during propagation
    end

    methods
        %% CONSTRUCTOR
        function obj = Sim(resolution, fieldLength, varargin)
            p = inputParser;
            addRequired(p, 'resolution', @isnumeric);
            addRequired(p, 'fieldLength', @isnumeric);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'lambda', 632e-9, @isnumeric);
            addParameter(p, 'paddingRatio', 1, @isnumeric);
            parse(p, resolution, fieldLength, varargin{:});

            % Assign parsed values to object properties
            obj.resolution = p.Results.resolution;
            obj.fieldLength = p.Results.fieldLength;
            obj.dim = p.Results.dim;
            obj.lambda = p.Results.lambda;
            obj.paddingRatio = p.Results.paddingRatio;

            % Set additional object properties
            obj.samples = round(obj.fieldLength / obj.resolution);
            obj.elements = {};
            obj.distances = [];
        end

        %% ELEMENT MANAGEMENT
        function addElement(obj, dist, element)
            % Add an optical element to the system with a specific distance
            % from the previous element.

            if nargin < 3
                dist = 0; % Default distance if not provided
            end

            if obj.dim ~= element.dim
                error('Dimensionality of element must match the system dimensionality.');
            end

            obj.elements{end+1} = element;
            obj.distances(end+1) = dist;
        end

        function addLens(obj, dist, focalLength, varargin)
            lens = Lens(focalLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, lens);
        end

        function addDiffuser(obj, dist, roughness, correlation_length, varargin)
            diffuser = Diffuser(roughness, correlation_length, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, diffuser);
        end

        function addAperture(obj, dist, varargin)
            aperture = Aperture('dim', obj.dim, varargin{:});
            obj.addElement(dist, aperture);
        end

        %% PROPAGATION METHODS
        % These methods are typically wrappers around the propagation calculation methods.
        function field = propToIndex(obj, field, targetIndex, varargin)
            % Propagate to a specific index of an element. Typically,
            % other propagation functions are wrappers of this function.
            % Note: The element at the index specified is applied to the
            % field.

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, varargin{:});
            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;

            % Propagate the field to a specific element by index
            if targetIndex > length(obj.elements)
                error('Target index exceeds the number of elements in the system.');
            end

            if verbose
                field.disp('Initial Field');
            end

            cumulativeDist = 0;

            % Loop through elements until target index is reached
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

                % Apply element's phase shift or aperture
                field = obj.elements{i}.apply(field);

                if verbose
                    fprintf('Propagated to element %d (%s) at distance: %.3e m\n', i, obj.elements{i}.name, cumulativeDist);
                    field.disp(sprintf('After element %d (%s)\nDist: %.3e m', i, obj.elements{i}.name, cumulativeDist));
                end
            end
        end

        function field = prop(obj, field, varargin)
            % Propagate through all elements in the system and optionally propagate a specified distance after the last element.

            % Parse input arguments
            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'Field'));
            addOptional(p, 'distAfter', 0, @isnumeric);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, varargin{:});

            % Extract parsed arguments
            distAfter = p.Results.distAfter;
            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;

            if verbose
                fprintf('Starting propagation through all elements in the system\n');
            end

            % Propagate through all elements
            lastElementIndex = length(obj.elements);
            if lastElementIndex > 0
                field = obj.propToIndex(field, lastElementIndex, 'verbose', verbose, 'propMethod', propMethod);
                currDist = sum(obj.distances);
            else
                currDist = 0;
            end

            % Propagate beyond the last element if specified
            if distAfter > 0
                if obj.dim == 1
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation1D(field, distAfter);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation1D(field, distAfter);
                    end
                else
                    if strcmp(propMethod, 'as')
                        field = obj.angularSpectrumPropagation(field, distAfter);
                    elseif strcmp(propMethod, 'rs')
                        field = obj.rayleighSommerfeldPropagation(field, distAfter);
                    end
                end
                currDist = currDist + distAfter;

                if verbose
                    fprintf('Propagated to distance: %.3e m\n', currDist);
                    field.disp(sprintf('Final Field\nDist: %.3e m', currDist));
                end
            end
        end

        function [field, currDist] = propToDist(obj, field, targetDist, varargin)
            % Propagate up to a specific distance in the system

            % Parse input arguments
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
                        field.disp('Field after propagation (no elements)');
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
                    field.disp('Field after propagation to target distance');
                end
            end
        end

        function field = propToElement(obj, field, targetName, varargin)
            % Propagate to an element by its name

            % Parse input arguments
            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'Field'));
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

            % Use the specified propagation method
            field = obj.propToIndex(field, elementIndex, 'verbose', verbose, 'propMethod', propMethod);

            if verbose
                cumulativeDist = sum(obj.distances(1:elementIndex));
                fprintf('Finished propagating to element named: %s at distance: %.3e m\n', targetName, cumulativeDist);
                field.disp(sprintf('Field after propagation to element (%s), Distance: %.3e m', targetName, cumulativeDist));
            end
        end

        %% PROPAGATION CALCULATION METHODS
        function field = angularSpectrumPropagation(obj, field, z)
            % Angular Spectrum Propagation (2D)
            n0 = 1;
            dx = obj.resolution;
            uin = field.getComplexField();

            uin_padded = obj.addPadding(uin);
            [Ny, Nx] = size(uin_padded);
            k = 2 * pi / obj.lambda;

            dfx = 1 / Nx / dx;
            fx = -Nx / 2 * dfx : dfx : (Nx / 2 - 1) * dfx;

            dfy = 1 / Ny / dx;
            fy = -Ny / 2 * dfy : dfy : (Ny / 2 - 1) * dfy;

            % Kernel for propagation
            p = fftshift(k * z * sqrt(n0^2 - obj.lambda^2 * (ones(Ny, 1) * (fx.^2) + (fy'.^2) * ones(1, Nx))));
            p = p - p(1, 1);
            kernel = exp(1i * p);

            ftu = kernel .* fft2(uin_padded);
            uout_padded = ifft2(ftu);

            uout = obj.removePadding(uout_padded);
            field = field.setComplexField(uout);
        end

        function field = angularSpectrumPropagation1D(obj, field, z)
            % Angular Spectrum Propagation (1D)
            dx = obj.resolution;
            uin = field.getComplexField();

            Nx = length(uin);
            k = 2 * pi / obj.lambda;

            dfx = 1 / (Nx * dx);
            fx = (-Nx / 2 : Nx / 2 - 1) * dfx;

            % Kernel for propagation in 1D
            kernel = exp(1i * k * z * sqrt(1 - (obj.lambda * fx).^2));
            ftu = fftshift(fft(uin));
            ftu = ftu .* kernel;
            uout = ifft(ifftshift(ftu));

            field = field.setComplexField(uout);
        end

        function field = rayleighSommerfeldPropagation(obj, field, z)
            % Rayleigh-Sommerfeld Propagation (2D)
            if z == 0
                field = field.setComplexField(field.getComplexField());
                return;
            end

            dx = obj.resolution;
            n0 = 1;
            k = 2 * pi / obj.lambda * n0;
            uin = field.getComplexField();

            % Add padding to the input field
            uin_padded = obj.addPadding(uin);
            [Ny, Nx] = size(uin_padded);

            % Generate the kernel grid
            xkernel = (-Nx / 2 : (Nx / 2 - 1)) * dx;
            ykernel = (-Ny / 2 : (Ny / 2 - 1)) * dx;
            [X, Y] = meshgrid(xkernel, ykernel);
            r = sqrt(X.^2 + Y.^2 + z^2);

            % Rayleigh-Sommerfeld Kernel for propagation
            if z > 0
                kernel = (2 * pi)^-1 * exp(1i * k * r) ./ r.^2 .* z .* (1 ./ r - 1i * k);
            else
                kernel = conj((2 * pi)^-1 * exp(1i * k * r) ./ r.^2 .* z .* (1 ./ r - 1i * k));
            end

            % Perform convolution in the spectrum domain
            uout_padded = ifftshift(ifft2(fft2(uin_padded) .* fft2(kernel))) * dx * dx;

            % Remove padding from the propagated field
            uout = obj.removePadding(uout_padded);
            field = field.setComplexField(uout);
        end

        function field = rayleighSommerfeldPropagation1D(obj, field, z)
            % Rayleigh-Sommerfeld Propagation (1D)
            if z == 0
                field = field.setComplexField(field.getComplexField());
                return;
            end

            dx = obj.resolution;
            n0 = 1; % Refractive index
            k = 2 * pi / obj.lambda * n0; % Wave number
            uin = field.getComplexField();

            % Pad the input field to match the kernel size
            uin_padded = obj.addPadding(uin);
            Nx = length(uin_padded);

            xkernel = (-Nx / 2 : Nx / 2 - 1) * dx;
            r = sqrt(xkernel.^2 + z^2);

            % Rayleigh-Sommerfeld Kernel for propagation in 1D
            if z > 0
                kernel = (2 * pi)^-1 * exp(1i * k * r) ./ r .* (1 - 1i * k * r);
            else
                kernel = conj((2 * pi)^-1 * exp(1i * k * r) ./ r .* (1 - 1i * k * r));
            end

            % Perform convolution in the spectrum domain
            uout_padded = ifft(fft(uin_padded) .* fft(kernel)) * dx * dx;

            % Crop the output to the original field size
            uout = obj.removePadding(uout_padded);

            % Set the propagated field
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
        function wrappedPhi = wrap(~, phi)
            % Wrap phase to the range [-pi, pi]
            wrappedPhi = atan2(sin(phi), cos(phi));
        end

        function field = newField(obj)
            % Create a new Field instance with system dimensions
            field = Field(obj.dim, obj.fieldLength, obj.resolution, obj.lambda);
        end
    end
end
