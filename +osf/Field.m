classdef Field
    properties
        amplitude    % Amplitude of the field (1D or 2D array)
        phase        % Phase of the field (1D or 2D array)
        fieldLength  % Length of the field in meters (scalar for 1D, [width, height] for 2D)
        resolution   % Spatial resolution in meters per pixel (scalar)
        size
        lambda       % Wavelength of light in meters
        dim          % Dimensionality of the field (1 or 2)

        cmap
    end

    methods
        function obj = Field(dim, fieldLength, resolution, lambda)
            if ~ismember(dim, [1, 2])
                error('Dimensionality must be either 1 or 2.');
            end

            obj.dim = dim;
            obj.fieldLength = fieldLength;
            obj.resolution = resolution;
            obj.lambda = lambda;
            obj.cmap = 'bone';

            samples = round(fieldLength / resolution);
            if dim == 1
                obj.amplitude = ones(1, samples);
                obj.phase = zeros(1, samples);
                obj.size = [1 samples]
            elseif dim == 2
                obj.amplitude = ones(samples, samples);
                obj.phase = zeros(samples, samples);
                obj.size = [samples samples]
            end
        end

        function complexField = getComplexField(obj)
            complexField = obj.amplitude .* exp(1i * obj.phase);
        end

        function obj = setComplexField(obj, field)
            obj.amplitude = abs(field);
            obj.phase = angle(field);
        end

        function obj = newField(obj)
            obj.amplitude = ones(obj.size);
            obj.phase = zeros(obj.size);
        end

        function obj = setAmplitude(obj, value, varargin)
            if nargin < 2
                error('setAmplitude requires at least a value argument.');
            end

            shape = 'global';
            region_size = [];
            position = [];

            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);

            if ~strcmp(type, 'global')
                if nargin < 4, error('Size is required for shape: %s', type); end
                region_size = varargin{2};
                if nargin > 4, position = varargin{3}; end
                region_size = round(region_size / obj.resolution);
                position = round(position / obj.resolution);
                mask = osf.utils.genMask(obj.size, type, region_size, position);
                obj.amplitude(mask) = value;
            else
                obj.amplitude(:) = value;
            end
        end

        function obj = addAmplitude(obj, value, varargin)
            if nargin < 2
                error('addAmplitude requires at least a value argument.');
            end

            shape = 'global';
            region_size = [];
            position = [];

            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);

            if ~strcmp(type, 'global')
                if nargin < 4, error('Size is required for shape: %s', type); end
                region_size = varargin{2};
                if nargin > 4, position = varargin{3}; end
                region_size = round(region_size / obj.resolution);
                position = round(position / obj.resolution);
                mask = osf.utils.genMask(obj.size, type, region_size, position);
                obj.amplitude(mask) = obj.amplitude(mask) + value;
            else
                obj.amplitude = obj.amplitude + value;
            end
        end

        function obj = setPhase(obj, value, varargin)
            if nargin < 2
                error('setPhase requires at least a value argument.');
            end

            shape = 'global';
            region_size = [];
            position = [];

            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);

            if ~strcmp(type, 'global')
                if nargin < 4, error('Size is required for shape: %s', type); end
                region_size = varargin{2};
                if nargin > 4, position = varargin{3}; end
                region_size = round(region_size / obj.resolution);
                position = round(position / obj.resolution);
                mask = osf.utils.genMask(obj.size, type, region_size, position);
                obj.phase(mask) = value;
            else
                obj.phase(:) = value;
            end
        end

        function obj = addPhase(obj, value, varargin)
            if nargin < 2
                error('addPhase requires at least a value argument.');
            end

            shape = 'global';
            region_size = [];
            position = [];

            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);

            if ~strcmp(type, 'global')
                if nargin < 4, error('Size is required for shape: %s', type); end
                region_size = varargin{2};
                if nargin > 4, position = varargin{3}; end
                region_size = round(region_size / obj.resolution);
                position = round(position / obj.resolution);
                mask = osf.utils.genMask(obj.size, type, region_size, position);
                obj.phase(mask) = obj.phase(mask) + value;
            else
                obj.phase = obj.phase + value;
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
                    pos = round(obj.size(1) / 2); % Middle row for x-axis cross section
                else
                    pos = round(obj.size(2) / 2); % Middle column for y-axis cross section
                end
            end

            % Extract cross-section based on axis type
            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, obj.size(2));
                amplitudeData = obj.amplitude(pos, :);
            elseif strcmp(axisType, 'y')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, obj.size(1));
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
                    pos = round(obj.size(1) / 2); % Middle row for x-axis cross section
                else
                    pos = round(obj.size(2) / 2); % Middle column for y-axis cross section
                end
            end

            % Extract cross-section based on axis type
            if strcmp(axisType, 'x')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, obj.size(2));
                phaseData = obj.phase(pos, :);
            elseif strcmp(axisType, 'y')
                xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, obj.size(1));
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
            targetSize = obj.size;
            img = imresize(img, targetSize, 'bilinear');

            % Assign to specified property
            if strcmp(type, 'phase')
                obj.phase = img;
            else
                obj.amplitude = img;
            end
        end

        function operationCheck(a, b)
            if ~isa(a, 'osf.Field') || ~isa(b, 'osf.Field')
                error('Both operands must be of class Field.');
            end
            if a.dim ~= b.dim
                error('Fields must have the same dimensionality.');
            end
            if ~isequal(a.size, b.size) || ~isequal(a.resolution, b.resolution)
                error('Fields must have the same size and resolution.');
            end
        end

        function result = plus(a, b)
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude + b.amplitude;
            result.phase = a.phase + b.phase;
        end

        function result = minus(a, b)
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude - b.amplitude;
            result.phase = a.phase - b.phase;
        end

        function result = times(a, b)
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude .* b.amplitude;
            result.phase = a.phase .* b.phase;
        end

        function result = mtimes(a, b)
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude .* b.amplitude;
            result.phase = a.phase .* b.phase;
        end

        function result = rdivide(a, b)
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude ./ b.amplitude;
            result.phase = a.phase ./ b.phase;
        end

    end
end
