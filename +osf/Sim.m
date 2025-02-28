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

            obj.resolution = p.Results.resolution;
            obj.fieldLength = p.Results.fieldLength;
            obj.dim = p.Results.dim;
            obj.lambda = p.Results.lambda;
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
            lens = osf.Lens(focalLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, lens);
        end

        function diffuser = addDiffuser(obj, dist, roughness, correlationLength, varargin)
            diffuser = osf.Diffuser(roughness, correlationLength, 'dim', obj.dim, varargin{:});
            obj.addElement(dist, diffuser);
        end

        function aperture = addAperture(obj, dist, varargin)
            aperture = osf.Aperture('dim', obj.dim, varargin{:});
            obj.addElement(dist, aperture);
        end

        function plane = addPlane(obj, dist, varargin)
            plane = osf.Plane('dim', obj.dim, varargin{:});
            obj.addElement(dist, plane);
        end

        %% PROPAGATION METHODS
        % These methods are typically wrappers around the propagation calculation methods.

        function field = propToIndex(obj, field, targetIndex, varargin)
            % Propagate to a specific index of an element. Typically,
            % other propagation methods are wrappers of this function.
            % Note: The element at the specified index is applied to the
            % field.

            p = inputParser;
            addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, varargin{:});
            verbose = p.Results.verbose;
            propMethod = p.Results.propMethod;

            if targetIndex > length(obj.elements)
                error('Target index exceeds the number of elements in the system.');
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

                % Apply element's phase shift and aperture
                field = obj.elements{i}.apply(field);

                if verbose
                    fprintf('Propagated to element %d (%s) at distance: %.3e m\n', i, obj.elements{i}.name, cumulativeDist);
                    field.disp('title', sprintf('After element %d (%s)\nDist: %.0f mm', i, obj.elements{i}.name, 1000*cumulativeDist));
                end
            end
        end

        function field = prop(obj, field, varargin)
            % Propagate through all elements in the system and optionally
            % propagate a specified distance after the last element.

            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'osf.Field'));
            addOptional(p, 'distAfter', 0, @isnumeric);
            addParameter(p, 'verbose', false, @islogical);
            addParameter(p, 'propMethod', 'as', @(x) ischar(x) && ismember(x, {'as', 'rs'}));
            parse(p, field, varargin{:});

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
                    field.disp('title', sprintf('Final Field\nDist: %.3e m', currDist));
                end
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
                        field.disp('title', 'Field after propagation (no elements)');
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
                    field.disp('title', 'Field after propagation to target distance');
                end
            end
        end

        function field = propToElement(obj, field, targetName, varargin)
            % Propagate to an element by its name

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

            field = obj.propToIndex(field, elementIndex, 'verbose', verbose, 'propMethod', propMethod);

            if verbose
                cumulativeDist = sum(obj.distances(1:elementIndex));
                fprintf('Finished propagating to element named: %s at distance: %.3e m\n', targetName, cumulativeDist);
                field.disp('title', sprintf('Field after propagation to element (%s), Distance: %.3e m', targetName, cumulativeDist));
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
            lambda = obj.lambda;

            % Frequency axis computation.
            fx = (-Nx/2:Nx/2-1) / (Nx * dx);
            fy = (-Ny/2:Ny/2-1) / (Ny * dy);
            [FX_sqr, FY_sqr] = meshgrid(fx.^2, fy.^2);

            % Compute longitudinal wavevector component k_z.
            k_z = 2 * pi * sqrt((1 / lambda^2) - FX_sqr - FY_sqr);
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
            k = 2 * pi / obj.lambda;

            dfx = 1 / (Nx * dx);
            fx = (-Nx / 2 : Nx / 2 - 1) * dfx;

            % Kernel for propagation in 1D
            kernel = exp(1i * k * z * sqrt(1 - (obj.lambda * fx).^2));
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

        function wrappedPhi = wrap(~, phi)
            % Wrap phase to the range [-pi, pi]
            wrappedPhi = atan2(sin(phi), cos(phi));
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

            field = osf.Field(obj.dim, obj.fieldLength, obj.resolution, obj.lambda);

            % Only override cmap if the user provides an input
            if ~isempty(p.Results.cmap)
                field.cmap = p.Results.cmap;
            end
        end

        function disp(obj)
            % Plots an overhead view of the optical setup, showing elements in order
            % with their respective distances.

            % Ensure there are elements to plot
            if isempty(obj.elements)
                error('No elements to display in the simulation.');
            end

            % Initialize figure with wider and slightly taller dimensions
            fig = figure('Color', 'white', 'Position', [600 300 1200 400]);
            hold on;
            axis equal;

            % Remove grid, tick marks, and box outline
            ax = gca;
            ax.XColor = 'none';
            ax.YColor = 'none';
            ax.XTick = [];
            ax.YTick = [];
            ax.Box = 'off';

            % Define common element height and label offset
            elementHeight = 0.04;  % Consistent height for all elements
            labelOffset = 0.015;    % Offset for labels

            % Initialize position tracker
            currentX = 0;
            componentCenters = []; % Store center positions for tick marks

            % Iterate over elements and update position first
            for i = 1:length(obj.elements)
                % Update position based on the corresponding distance before plotting the element
                currentX = currentX + obj.distances(i);
                element = obj.elements{i};

                % Plot the element based on its type
                switch lower(element.elementType)
                    case 'lens'
                        plotLens(currentX);
                    case 'diffuser'
                        plotDiffuser(currentX);
                    case 'plane'
                        plotPlane(currentX);
                    otherwise
                        warning('Unknown element type: %s', element.elementType);
                        plotUnknown(currentX);
                end

                % Add a label above the element
                text(currentX+.005, elementHeight/2 + labelOffset, element.name, ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, ...
                    'FontWeight', 'bold', 'Rotation', 45);

                % Store component center position
                componentCenters = [componentCenters, currentX];
            end

            % Adjust axis limits (unchanged settings)
            xlim([min(0, -.5*currentX), currentX + .5*currentX]);
            ylim([-0.1, 0.1]);

            % Draw a thick black baseline and tick marks with distance labels
            baselineY = -0.03;
            plot([min(componentCenters), max(componentCenters)], [baselineY, baselineY], 'k', 'LineWidth', 2);
            tickHeight = 0.002;
            for i = 1:length(componentCenters)
                xTick = componentCenters(i);
                plot([xTick, xTick], [baselineY, baselineY + tickHeight], 'k', 'LineWidth', 2);
                if i < length(componentCenters)
                    midX = (componentCenters(i) + componentCenters(i+1)) / 2;
                    text(midX, baselineY - 0.01, sprintf('%.0fmm', obj.distances(i+1)*1000), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end

            % --- Paraxial Ray Overlay ---
            % Create a ParaxialSystem object from this Sim object and solve it.
            parax = osf.parax.ParaxialSystem(obj);
            parax = osf.parax.ParaxialSystem.solveSystem(parax);

            % Compute the cumulative distances along the system.
            dist = cumsum(parax.distances);

            % Overlay the computed marginal and chief rays:
            plot(dist, parax.marginalRay.heights*(0.8*elementHeight/max(2*abs(parax.marginalRay.heights))), 'r', 'LineWidth', 0.5);
            plot(dist, -parax.marginalRay.heights*(0.8*elementHeight/max(2*abs(parax.marginalRay.heights))), 'r', 'LineWidth', 0.5);
            plot(dist, parax.chiefRay.heights*(0.8*elementHeight/max(2*abs(parax.chiefRay.heights))), 'b', 'LineWidth', 0.5);
            plot(dist, -parax.chiefRay.heights*(0.8*elementHeight/max(2*abs(parax.chiefRay.heights))), 'b', 'LineWidth', 0.5);

            hold off;

            % --- Nested plotting functions ---
            function plotLens(x)
                width = 0.005;
                t = linspace(0, pi, 20);
                X = [cos(t), -cos(t)] * width/2 + x;
                Y = [sin(t), -sin(t)] * elementHeight/2;
                fill(X, Y, [.8 .8 .8], 'EdgeColor', 'k', 'LineWidth', 2);
            end

            function plotDiffuser(x)
                width = 0.002;
                rectangle('Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
                    'FaceColor', 'g', 'EdgeColor', 'k', 'LineWidth', 1.5);
            end

            function plotPlane(x)
                plot([x, x], [-elementHeight/2, elementHeight/2], 'k', 'LineWidth', 2);
            end

            function plotUnknown(x)
                width = 0.003;
                rectangle('Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
                    'FaceColor', 'r', 'EdgeColor', 'k', 'LineWidth', 1.5);
            end

            % To make sure the disp shows before some other process starts running.
            pause(.001);
        end

        function print(obj)
            % PRINT Prints simulation parameters in a nicely formatted manner.
            %   Length measurements are converted to mm (or um for resolution)
            %   with no digits after the decimal, and element names are printed in a 
            %   bracketed list. If an element's name is empty, its type is used instead.

            distances_mm   = obj.distances * 1e3;
            fieldLength_mm = obj.fieldLength * 1e3;
            resolution_um  = obj.resolution * 1e6;
            lambda_nm      = obj.lambda * 1e9;

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
            fprintf('  Wavelength:     %d nm\n', lambda_nm);
            fprintf('  Padding Ratio:  %.1f\n\n', obj.paddingRatio);
        end

    end
end
