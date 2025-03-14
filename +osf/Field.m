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

        function W = wdf(obj)

            field.disp(cross=true);
            [xdata, crossdata] = field.getPhaseCross();
            [d, f, t] = wvd(crossdata, 1000000, 'smoothedPseudo', MinThreshold=0, NumFrequencyPoints=1000);
            figure; imagesc(t,f,d(1:100,:)); ylabel('Local Spatial Frequency (m^{-1})'); xlabel('x (m)');

        end

        function obj = fft(obj)
            obj.setComplexField(fftshift(fft(obj.getComplexField())));
        end

        function [xAxis, amplitudeData] = getAmplitudeCross(obj, axisType, pos)
            % getAmplitudeCross Returns the amplitude cross section of a 2D field.
            % Default axis is 'x', and default position is the middle of the field.
            % Throws an error if the field is not 2D.

            if obj.dim ~= 2
                error('getAmplitudeCross is defined only for 2D fields.');
            end

            % Set default values if not provided
            if nargin < 2 || isempty(axisType)
                axisType = 'x'; % Default to x-axis
            end
            if nargin < 3 || isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(size(obj.amplitude, 1) / 2); % Middle row for x-axis cross section
                else
                    pos = round(size(obj.amplitude, 2) / 2); % Middle column for y-axis cross section
                end
            end

            % Extract cross-section based on axis type
            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.amplitude, 2));
                amplitudeData = obj.amplitude(pos, :);
            elseif strcmp(axisType, 'y')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.amplitude, 1));
                amplitudeData = obj.amplitude(:, pos);
            else
                error('Invalid axisType. Must be ''x'' or ''y''.');
            end
        end

        function [xAxis, phaseData] = getPhaseCross(obj, axisType, pos)
            % getPhaseCross Returns the phase cross section of a 2D field.
            % Default axis is 'x', and default position is the middle of the field.
            % Throws an error if the field is not 2D.

            if obj.dim ~= 2
                error('getPhaseCross is defined only for 2D fields.');
            end

            % Set default values if not provided
            if nargin < 2 || isempty(axisType)
                axisType = 'x'; % Default to x-axis
            end
            if nargin < 3 || isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(size(obj.phase, 1) / 2); % Middle row for x-axis cross section
                else
                    pos = round(size(obj.phase, 2) / 2); % Middle column for y-axis cross section
                end
            end

            % Extract cross-section based on axis type
            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.phase, 2));
                phaseData = obj.phase(pos, :);
            elseif strcmp(axisType, 'y')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, size(obj.phase, 1));
                phaseData = obj.phase(:, pos);
            else
                error('Invalid axisType. Must be ''x'' or ''y''.');
            end
        end

        function obj = unwrap(obj)
            obj.phase = osf.utils.phaseUnwrap(obj.phase);
        end

        function obj = show(obj, varargin)
            p = inputParser;
            addParameter(p, 'unwrap', true, @islogical);
            p.KeepUnmatched = true;

            parse(p, varargin{:});
            unwrapFlag = p.Results.unwrap;

            phaseType = 'phase';
            if unwrapFlag
                phaseType = 'phase.unwrap';
            end

            osf.Dashboard([2,1], p.Unmatched).show(obj, 'amplitude').show(obj, phaseType);
        end

        function obj = img(obj, path, varargin)
            p = inputParser;
            addOptional(p, 'type', 'phase', @(x) ismember(x, {'phase', 'amplitude'}));
            addOptional(p, 'range', [0 1], @(x) isnumeric(x) && numel(x) == 2);
            parse(p, varargin{:});

            type = p.Results.type;
            range = p.Results.range;

            % Read image and convert to grayscale if necessary
            img = imread(path);
            if size(img, 3) == 3
                img = rgb2gray(img); % Convert to grayscale if RGB
            end
            img = double(img); % Convert to double for rangealization

            % Normalize image to the specified range
            img = (img - min(img(:))) / (max(img(:)) - min(img(:))); % Scale to [0,1]
            img = range(1) + img * (range(2) - range(1)); % Scale to [range(1), range(2)]

            % Resize image to match field dimensions
            targetSize = size(obj.amplitude);
            img = imresize(img, targetSize, 'bilinear');

            % Assign to specified property
            if strcmp(type, 'phase')
                obj.phase = img;
            else
                obj.amplitude = img;
            end
        end

    end
end
