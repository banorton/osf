classdef Field
    properties
        amplitude    % Amplitude of the field (1D or 2D array)
        phase        % Phase of the field (1D or 2D array)
        fieldLength  % Length of the field in meters (scalar for 1D, [width, height] for 2D)
        resolution   % Spatial resolution in meters per pixel (scalar)
        lambda       % Wavelength of light in meters
        dim          % Dimensionality of the field (1 or 2)

        cmap
    end

    methods
        function obj = Field(dim, fieldLength, resolution, lambda)
            % Validate dimensionality
            if ~ismember(dim, [1, 2])
                error('Dimensionality must be either 1 or 2.');
            end

            obj.dim = dim;
            obj.fieldLength = fieldLength;
            obj.resolution = resolution;
            obj.lambda = lambda;
            obj.cmap = 'bone';

            % Initialize amplitude and phase based on dimensionality
            if dim == 1
                samples = round(fieldLength / resolution);
                obj.amplitude = zeros(1, samples);
                obj.phase = zeros(1, samples);
            elseif dim == 2
                samples = round(fieldLength / resolution);
                obj.amplitude = zeros(samples, samples);
                obj.phase = zeros(samples, samples);
            end
        end

        function complexField = getComplexField(obj)
            % Returns the complex representation of the field
            complexField = obj.amplitude .* exp(1i * obj.phase);
        end

        function obj = setComplexField(obj, field)
            % Set the amplitude and phase properties from a complex field
            obj.amplitude = abs(field);
            obj.phase = angle(field);
        end

        function obj = resetField(obj)
            % Reset the field amplitude to all 1s and phase to all 0s
            obj.amplitude = ones(size(obj.amplitude));
            obj.phase = zeros(size(obj.phase));
        end

        function obj = addPhase(obj, value, varargin)
            % addPhase applies a phase shift to the field.
            %
            % Usage:
            %   field.addPhase(value)                         % Global phase shift
            %   field.addPhase(value, 'rectangle', size)     % Rectangular shift
            %   field.addPhase(value, 'circle', radius)      % Circular shift
            %   field.addPhase(value, 'annulus', [r1, r2])   % Annular shift
            %   field.addPhase(value, 'rectangle', size, offset) % Rect at offset
            %   field.addPhase(value, 'circle', radius, offset)  % Circle at offset
            %   field.addPhase(value, 'annulus', [r1, r2], offset) % Annulus at offset

            % Ensure first argument (value) is provided
            if nargin < 2
                error('addPhase requires at least a value argument.');
            end

            % Default values
            shape = 'global';  % Default shape type
            region_size = [];
            position = [];

            % Parse the second argument (shape/type)
            if nargin > 2
                shape = varargin{1};
            end

            % Standardize shape type
            type = obj.standardizeType(shape);

            % If 'global', no additional arguments are allowed
            if strcmp(type, 'global')
                if nargin > 3
                    error('Global phase shift does not take a size or position argument.');
                end
            else
                % If type is not 'global', the 3rd argument (size) is required
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end

                % Assign and validate the size argument
                region_size = varargin{2};
                if strcmp(type, 'circle') && ~isscalar(region_size)
                    error('Circle size must be a scalar (radius).');
                end
                if strcmp(type, 'rectangle') && isscalar(region_size) && obj.dim == 2
                    region_size = [region_size, region_size]; % Convert scalar to [w, h] if 2D
                end
                if strcmp(type, 'annulus') && (~isnumeric(region_size) || numel(region_size) ~= 2)
                    error('Annulus size must be a two-element vector [inner_radius, outer_radius].');
                end
            end

            % Assign position offset if provided
            if nargin > 4
                position = varargin{3};
            end

            region_size = round(region_size / obj.resolution);
            position = round(position / obj.resolution);

            % Apply modification (dimensionality checks happen inside applyMod)
            if obj.dim == 1
                obj.phase = obj.applyModification1D(obj.phase, value, type, region_size, position);
            elseif obj.dim == 2
                obj.phase = obj.applyModification2D(obj.phase, value, type, region_size, position);
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function obj = addAmplitude(obj, value, varargin)
            % addAmplitude applies an amplitude change to the field.
            %
            % Usage:
            %   field.addAmplitude(value)                         % Global amplitude change
            %   field.addAmplitude(value, 'rectangle', size)     % Rectangular shift
            %   field.addAmplitude(value, 'circle', radius)      % Circular shift
            %   field.addAmplitude(value, 'annulus', [r1, r2])   % Annular shift
            %   field.addAmplitude(value, 'rectangle', size, offset) % Rect at offset
            %   field.addAmplitude(value, 'circle', radius, offset)  % Circle at offset
            %   field.addAmplitude(value, 'annulus', [r1, r2], offset) % Annulus at offset

            % Ensure first argument (value) is provided
            if nargin < 2
                error('addAmplitude requires at least a value argument.');
            end

            % Default values
            shape = 'global';  % Default shape type
            region_size = [];
            position = [];

            % Parse the second argument (shape/type)
            if nargin > 2
                shape = varargin{1};
            end

            % Standardize shape type
            type = obj.standardizeType(shape);

            % If 'global', no additional arguments are allowed
            if strcmp(type, 'global')
                if nargin > 3
                    error('Global amplitude shift does not take a size or position argument.');
                end
            else
                % If type is not 'global', the 3rd argument (size) is required
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end

                % Assign and validate the size argument
                region_size = varargin{2};
                if strcmp(type, 'circle') && ~isscalar(region_size)
                    error('Circle size must be a scalar (radius).');
                end
                if strcmp(type, 'rectangle') && isscalar(region_size) && obj.dim == 2
                    region_size = [region_size, region_size]; % Convert scalar to [w, h] if 2D
                end
                if strcmp(type, 'annulus') && (~isnumeric(region_size) || numel(region_size) ~= 2)
                    error('Annulus size must be a two-element vector [inner_radius, outer_radius].');
                end
            end

            % Assign position offset if provided
            if nargin > 4
                position = varargin{3};
            end

            region_size = round(region_size / obj.resolution);
            position = round(position / obj.resolution);

            % Apply modification (dimensionality checks happen inside applyMod)
            if obj.dim == 1
                obj.amplitude = obj.applyModification1D(obj.amplitude, value, type, region_size, position);
            elseif obj.dim == 2
                obj.amplitude = obj.applyModification2D(obj.amplitude, value, type, region_size, position);
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function type = standardizeType(~, inputType)
            % Maps multiple aliases to a single canonical shape name.
            
            % Define mapping of all possible aliases
            typeMap = struct(...
                'global', 'global', 'g', 'global', 'glob', 'global', ...
                'rect', 'rectangle', 'rectangle', 'rectangle', 'r', 'rectangle', 'rectangular', 'rectangle', ...
                'circ', 'circle', 'circle', 'circle', 'c', 'circle', 'circular', 'circle', ...
                'annu', 'annulus', 'annulus', 'annulus', 'a', 'annulus', 'annular', 'annulus');
        
            % Ensure input is a string and convert to lowercase
            if ~ischar(inputType) && ~isstring(inputType)
                error('Input type must be a string.');
            end
            inputType = lower(char(inputType)); % Convert string to char if needed
        
            % Check if input exists in typeMap
            if isfield(typeMap, inputType)
                type = typeMap.(inputType);
            else
                error('Invalid type: %s. Must be ''global'', ''rectangle'', ''circle'', or ''annulus''.', inputType);
            end
        end

        function field = applyModification1D(~, field, value, type, region_size, position)
            % Applies a phase or amplitude change to a 1D field based on shape type.
            if isempty(position)
                position = 0;
            end

            len = length(field);
            middle = round(len / 2);  % Center of the field
            pos_shift = round(position / 2);  % Adjusted for pixel index shift

            switch type
            case 'global'
                field = field + value;

            case 'rectangle'
                if isempty(region_size)
                    error('Rectangle size must be specified for 1D.');
                end
                rect_half = round(region_size / 2);
                start_idx = max(1, middle + pos_shift - rect_half);
                end_idx = min(len, middle + pos_shift + rect_half);
                field(start_idx:end_idx) = field(start_idx:end_idx) + value;

            case 'circle'
                error('Circle modifications are not valid for 1D fields.');

            case 'annulus'
                error('Annulus modifications are not valid for 1D fields.');

            otherwise
                error('Unsupported type for 1D modification: %s', type);
            end
        end

        function field = applyModification2D(~, field, value, type, region_size, position)
            % Applies a phase or amplitude change to a 2D field based on shape type.
            if isempty(position)
                position = [0 0];
            end

            [height, width] = size(field);
            middle_x = round(width / 2);
            middle_y = round(height / 2);

            % Apply position shift
            pos_shift_x = round(position(1) / 2);
            pos_shift_y = round(position(2) / 2);

            center_x = middle_x + pos_shift_x;
            center_y = middle_y + pos_shift_y;

            switch type
            case 'global'
                field = field + value;

            case 'rectangle'
                if isempty(region_size)
                    error('Rectangle size must be specified.');
                end
                if isscalar(region_size)
                    region_size = [region_size, region_size]; % Convert scalar to [w, h]
                end
                half_w = round(region_size(1) / 2);
                half_h = round(region_size(2) / 2);
                start_x = max(1, center_x - half_w);
                end_x = min(width, center_x + half_w);
                start_y = max(1, center_y - half_h);
                end_y = min(height, center_y + half_h);
                field(start_y:end_y, start_x:end_x) = field(start_y:end_y, start_x:end_x) + value;

            case 'circle'
                if isempty(region_size) || ~isscalar(region_size)
                    error('Circle size must be a scalar (radius).');
                end
                radius = round(region_size);
                [X, Y] = meshgrid(1:width, 1:height);
                mask = ((X - center_x).^2 + (Y - center_y).^2) <= radius^2;
                field(mask) = field(mask) + value;

            case 'annulus'
                if isempty(region_size) || numel(region_size) ~= 2
                    error('Annulus size must be a two-element vector [inner_radius, outer_radius].');
                end
                outer_radius = round(region_size(1));
                inner_radius = round(region_size(2));
                [X, Y] = meshgrid(1:width, 1:height);
                mask = ((X - center_x).^2 + (Y - center_y).^2) <= outer_radius^2 & ...
                ((X - center_x).^2 + (Y - center_y).^2) > inner_radius^2;
                field(mask) = field(mask) + value;

            otherwise
                error('Unsupported type for 2D modification: %s', type);
            end
        end

        function wdfValues = wdf(obj)
            % Calculate the gradient or difference of the field phase
            if obj.dim == 2
                % 2D case: compute gradient of the phase
                [dx, dy] = gradient(obj.phase, obj.resolution);
                wdfValues = sqrt(dx.^2 + dy.^2); % Magnitude of the gradient
            elseif obj.dim == 1
                % 1D case: compute difference of the phase
                wdfValues = diff(obj.phase) / obj.resolution;
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function fftField = fft(obj)
            % Calculate the Fourier tranform of the complex field
            complexField = obj.getComplexField();
            if obj.dim == 1
                fftField = fftshift(fft(complexField));
            elseif obj.dim == 2
                fftField = fftshift(fft2(complexField));
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function [xAxis, amplitudeData] = getAmplitudeCross(obj, axisType, pos)
            % getAmplitudeCross Returns the amplitude cross section of a 2D field.
            % Throws an error if the field is not 2D.
            if obj.dim ~= 2
                error('getAmplitudeCross is defined only for 2D fields.');
            end

            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.amplitude,2));
                amplitudeData = obj.amplitude(pos, :);
            else  % axisType == 'y'
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.amplitude,1));
                amplitudeData = obj.amplitude(:, pos);
            end
        end

        function [xAxis, phaseData] = getPhaseCross(obj, axisType, pos)
            % getPhaseCross Returns the phase cross section of a 2D field.
            % Throws an error if the field is not 2D.
            if obj.dim ~= 2
                error('getPhaseCross is defined only for 2D fields.');
            end

            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.phase,2));
                phaseData = obj.phase(pos, :);
            else  % axisType == 'y'
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.phase,1));
                phaseData = obj.phase(:, pos);
            end
        end

        function fig = disp(obj, varargin)
            % DISP Displays the field amplitude and phase.
            %
            % By default (cross = false), it shows full-field amplitude and phase in a 2-row layout.
            % If 'cross' is true, it creates a 2x2 grid with full-field images in the left column and
            % cross sections in the right column.
            %
            % Optional name/value pairs:
            %   'titleName' - Overall figure title (default: "")
            %   'cross'     - Logical flag; if true, display additional cross sections (default: false)
            %   'axis'      - 'x' (default) or 'y' for the cross-section direction
            %   'pos'       - Row or column index for the cross section (default: center)

            p = inputParser;
            addParameter(p, 'titleName', "", @ischar);
            addParameter(p, 'cross', false, @islogical);
            addParameter(p, 'axis', 'x', @(x) ischar(x) && ismember(lower(x), {'x','y'}));
            addParameter(p, 'pos', [], @(x) isnumeric(x) && isscalar(x));
            parse(p, varargin{:});

            titleName = p.Results.titleName;
            crossFlag = p.Results.cross;
            axisType = lower(p.Results.axis);
            pos = p.Results.pos;

            % This function is defined only for 2D fields.
            if obj.dim ~= 2
                error('disp with cross option is defined only for 2D fields.');
            end

            [rows, cols] = size(obj.amplitude);
            if isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(rows/2);
                else
                    pos = round(cols/2);
                end
            end

            if crossFlag
                % Create a new figure with a 2x2 grid.
                fig = figure('Position', [300 100 1200 700]);

                % Left column: full-field images.
                xAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols)*1e3;
                yAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows)*1e3;

                % Top left: full amplitude image.
                subplot(2,2,1);
                imagesc(xAxis_full, yAxis_full, obj.amplitude);
                colormap(gca, 'gray');
                colorbar;
                title('Field Amplitude');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                % Bottom left: full phase image.
                subplot(2,2,3);
                imagesc(xAxis_full, yAxis_full, osf.utils.phase_unwrap(obj.phase));
                colormap(gca, obj.cmap);
                colorbar;
                title('Field Phase');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                % Right column: cross sections.
                % Amplitude cross.
                [xAxis_amp, ampCrossData] = getAmplitudeCross(obj, axisType, pos);
                subplot(2,2,2);
                plot(xAxis_amp*1e3, ampCrossData, 'LineWidth', 1.5);
                title(sprintf('Amplitude Cross (%s-axis, pos = %.2fmm)', axisType, pos*obj.resolution*1e3));
                xlabel('Position (mm)');
                ylabel('Amplitude');
                grid off;

                % Phase cross.
                [xAxis_phase, phaseCrossData] = getPhaseCross(obj, axisType, pos);
                subplot(2,2,4);
                plot(xAxis_phase*1e3, unwrap(phaseCrossData), 'LineWidth', 1.5);
                title(sprintf('Phase Cross (%s-axis, pos = %.2fmm)', axisType, pos*obj.resolution*1e3));
                xlabel('Position (mm)');
                ylabel('Phase (rad)');
                grid off;

            else
                % Without cross sections: full-field images in a 2-row layout.
                fig = figure('Position', [794 145 391 725]);
                xAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols)*1e3;
                yAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows)*1e3;

                subplot(2,1,1);
                imagesc(xAxis_full, yAxis_full, obj.amplitude);
                colormap(gca, 'gray');
                colorbar;
                title('Field Amplitude');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                subplot(2,1,2);
                imagesc(xAxis_full, yAxis_full, osf.utils.phase_unwrap(obj.phase));
                colorbar;
                colormap(gca, obj.cmap);
                title('Field Phase');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;
            end

            if ~isempty(titleName)
                st = sgtitle(titleName, 'FontSize', 18, 'FontWeight', 'bold');
                st.Tag = 'text';
            end

            obj.applyTheme(fig);

            pause(.001);
        end

        function fig = cross(obj, varargin)
            % CROSS Displays the amplitude and phase of a 2D field, with optional cross sections.
            %
            % This function is defined only for 2D fields.
            % It shows full-field amplitude and phase images, and, if the 'cross' flag
            % is true, it also displays the amplitude and phase cross sections in a 2x2 grid.
            %
            % Optional name/value pairs:
            %   'display' - 'both' (default), 'amplitude', or 'phase'
            %   'axis'    - 'x' (default) or 'y' for the cross-section direction
            %
            %   'pos'     - Row or column index for the cross section (default: center)
            %   'cross'   - Logical flag; if true, display the cross sections (default: false)

            p = inputParser;
            addParameter(p, 'display', 'both', @(x) ischar(x) && ismember(lower(x), {'both','amplitude','phase'}));
            addParameter(p, 'axis', 'x', @(x) ischar(x) && ismember(lower(x), {'x','y'}));
            addParameter(p, 'pos', [], @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'cross', false, @islogical);
            parse(p, varargin{:});

            displayChoice = lower(p.Results.display);
            axisType = lower(p.Results.axis);
            pos = p.Results.pos;
            crossFlag = p.Results.cross;

            % Set default pos value if not provided.
            [rows, cols] = size(obj.amplitude);
            if isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(rows/2);
                else
                    pos = round(cols/2);
                end
            end

            if crossFlag
                % Create a new figure with a 2x2 grid.
                fig = figure('Position', [300 100 1200 700]);

                % Left column: full-field images.
                xAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols)*1e3;
                yAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows)*1e3;

                % Top left: full-field amplitude.
                subplot(2,2,1);
                imagesc(xAxis_full, yAxis_full, obj.amplitude);
                colormap(gca, 'gray');
                colorbar;
                title('Field Amplitude');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                % Bottom left: full-field phase.
                subplot(2,2,3);
                imagesc(xAxis_full, yAxis_full, osf.utils.phase_unwrap(obj.phase));
                colormap(gca, obj.cmap);
                colorbar;
                title('Field Phase');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                % Right column: cross sections using external helper functions.
                % Amplitude cross section.
                [xAxis_amp, ampCrossData] = getAmplitudeCross(obj, axisType, pos);
                subplot(2,2,2);
                plot(xAxis_amp*1e3, ampCrossData, 'LineWidth', 1.5);
                title(sprintf('Amplitude Cross (%s-axis, pos = %.2fmm)', axisType, pos*obj.resolution*1e3));
                xlabel('Position (mm)');
                ylabel('Amplitude');
                grid off;

                % Phase cross section.
                [xAxis_phase, phaseCrossData] = getPhaseCross(obj, axisType, pos);
                subplot(2,2,4);
                plot(xAxis_phase*1e3, unwrap(phaseCrossData), 'LineWidth', 1.5);
                title(sprintf('Phase Cross (%s-axis, pos = %.2fmm)', axisType, pos*obj.resolution*1e3));
                xlabel('Position (mm)');
                ylabel('Phase (rad)');
                grid off;
            else
                % Without cross sections: use a 2-row layout with full-field images.
                fig = figure('Position', [794 145 391 725]);
                xAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols);
                yAxis_full = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows);

                subplot(2,1,1);
                imagesc(xAxis_full, yAxis_full, obj.amplitude);
                colormap(gca, 'gray');
                colorbar;
                title('Field Amplitude');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;

                subplot(2,1,2);
                imagesc(xAxis_full, yAxis_full, osf.utils.phase_unwrap(obj.phase));
                colorbar;
                colormap(gca, obj.cmap);
                title('Field Phase');
                xlabel('x (mm)');
                ylabel('y (mm)');
                axis equal; axis tight;
            end

            obj.applyTheme(fig);

            pause(.001);
        end

        function dispWDF(obj)
            % Displays the WDF (gradient or difference) of the field phase
            wdfValues = obj.wdf();

            if obj.dim == 1
                samples = length(wdfValues);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                figure('Color', 'white', 'Position', [657 76 486 846]);
                plot(xAxis, wdfValues, 'LineWidth', 1.5);
                title('WDF of Field Phase');
                xlabel('Position (m)');
                ylabel('WDF Value');
                grid on;

            elseif obj.dim == 2
                samples = size(wdfValues, 1);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);
                yAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                figure('Color', 'white', 'Position', [657 76 486 846]);
                imagesc(xAxis, yAxis, wdfValues);
                colormap(gca, 'default');
                colorbar;
                title('WDF of Field Phase');
                xlabel('x (m)');
                ylabel('y (m)');
                axis equal;
                axis tight;

            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function dispFFT(obj)
            % Display the FFT of the field
            fftField = obj.fft();

            if obj.dim == 1
                samples = length(fftField);
                spatialFreq = linspace(-1 / (2 * obj.resolution), 1 / (2 * obj.resolution), samples);

                figure('Color', 'black', 'Position', [657 76 486 846]);

                subplot(2, 1, 1);
                plot(spatialFreq, abs(fftField), 'LineWidth', 1.5);
                title('FFT Amplitude');
                xlabel('Spatial Frequency (1/m)');
                ylabel('Amplitude');
                grid on;

                subplot(2, 1, 2);
                plot(spatialFreq, angle(fftField), 'LineWidth', 1.5);
                title('FFT Phase');
                xlabel('Spatial Frequency (1/m)');
                ylabel('Phase (radians)');
                grid on;

            elseif obj.dim == 2
                samples = size(fftField, 1);
                spatialFreq = linspace(-1 / (2 * obj.resolution), 1 / (2 * obj.resolution), samples);

                figure('Color', 'black', 'Position', [657 76 486 846]);

                subplot(2, 1, 1);
                imagesc(spatialFreq, spatialFreq, abs(fftField));
                colormap(gca, 'gray');
                colorbar;
                title('FFT Amplitude');
                xlabel('Spatial Frequency x (1/m)');
                ylabel('Spatial Frequency y (1/m)');
                axis equal;
                axis tight;

                subplot(2, 1, 2);
                imagesc(spatialFreq, spatialFreq, angle(fftField));
                colormap(gca, 'default');
                colorbar;
                title('FFT Phase');
                xlabel('Spatial Frequency x (1/m)');
                ylabel('Spatial Frequency y (1/m)');
                axis equal;
                axis tight;

            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function applyTheme(obj, fig)
            % Applies a consistent theme to a figure with a black background,
            % white text, and appropriate formatting.

            if nargin < 2
                fig = gcf; % Use current figure if none is specified
            end

            % Set figure background color
            fig.Color = 'black';

            % Find all axes in the figure and apply theme
            ax = findall(fig, 'Type', 'axes');
            for i = 1:length(ax)
                ax(i).XColor = 'white';  % Change X-axis color
                ax(i).YColor = 'white';  % Change Y-axis color
                ax(i).ZColor = 'white';  % Change Z-axis color (if 3D)
                ax(i).Color = 'black';   % Set axes background to black
                ax(i).GridColor = [0.5, 0.5, 0.5]; % Light gray grid lines
                ax(i).LineWidth = 1.2; % Make axis lines slightly thicker
            end

            % Change title, labels, and other text objects
            textObjs = findall(fig, 'Type', 'text');
            for i = 1:length(textObjs)
                textObjs(i).Color = 'white';  % Change text color
                textObjs(i).FontWeight = 'bold';
            end

            % Also set sgtitle objects (which have the Tag 'text') to white.
            sgObjs = findall(fig, 'Tag', 'text');
            for i = 1:length(sgObjs)
                sgObjs(i).Color = 'white';
            end

            % Apply to colorbars if present
            colorbars = findall(fig, 'Type', 'colorbar');
            for i = 1:length(colorbars)
                colorbars(i).Color = 'white';
            end
        end

    end
end
