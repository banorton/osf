classdef Field
    % Field class for representing an optical field.
    % Stores amplitude, phase, and spatial properties, and provides methods
    % to manipulate and visualize the field.
    
    properties
        amplitude    % Amplitude of the field (1D or 2D array)
        phase        % Phase of the field (1D or 2D array)
        fieldLength  % Field length in meters (scalar for 1D, [width, height] for 2D)
        resolution   % Spatial resolution in meters per pixel (scalar)
        size         % Size of the field matrix
        lambda       % Wavelength of light in meters
        dim          % Dimensionality of the field (1 or 2)
        cmap         % Colormap for visualization
    end
    
    properties (Dependent)
        intensity    % Intensity computed as |complexField|^2
        complexField % Complex field computed as amplitude .* exp(1i*phase)
    end
    
    methods
        function obj = Field(dim, fieldLength, resolution, lambda)
            % Constructs a new Field object.
            %
            % Required
            %   dim : Dimensionality (1 or 2)
            %   fieldLength : Field length in meters (scalar for 1D, [width, height] for 2D)
            %   resolution : Spatial resolution in meters per pixel
            %   lambda : Wavelength in meters
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
                obj.size = [1 samples];
            else
                obj.amplitude = ones(samples, samples);
                obj.phase = zeros(samples, samples);
                obj.size = [samples samples];
            end
        end
        
        function complexField = getComplexField(obj)
            % Returns the complex field computed as amplitude .* exp(1i*phase).
            %
            % Required
            %   obj : Field object
            complexField = obj.amplitude .* exp(1i * obj.phase);
        end
        
        function obj = setComplexField(obj, field)
            % Sets the Field object's amplitude and phase from the input complex field.
            %
            % Required
            %   field : A complex field matrix.
            obj.amplitude = abs(field);
            obj.phase = angle(field);
        end
        
        function obj = newField(obj)
            % Resets the Field object to default values (amplitude ones, phase zeros).
            %
            % Required
            %   obj : Field object
            obj.amplitude = ones(obj.size);
            obj.phase = zeros(obj.size);
        end
        
        function obj = setAmplitude(obj, value, varargin)
            % Sets the amplitude of the field.
            %
            % Required
            %   value : New amplitude value.
            %
            % Optionals
            %   shape : Shape specification (default 'global').
            %
            % Parameters
            %   region_size : Size of the region (if shape is not 'global').
            %   position : Position offset (optional).
            if nargin < 2
                error('setAmplitude requires at least a value argument.');
            end
            shape = 'global';
            region_size = [];
            position = [];
            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);
            if ~strcmp(type, 'global')
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end
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
            % Adds a value to the current amplitude of the field.
            %
            % Required
            %   value : Value to add.
            %
            % Optionals
            %   shape : Shape specification (default 'global').
            %
            % Parameters
            %   region_size : Size of the region (if shape is not 'global').
            %   position : Position offset (optional).
            if nargin < 2
                error('addAmplitude requires at least a value argument.');
            end
            shape = 'global';
            region_size = [];
            position = [];
            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);
            if ~strcmp(type, 'global')
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end
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
            % Sets the phase of the field.
            %
            % Required
            %   value : New phase value.
            %
            % Optionals
            %   shape : Shape specification (default 'global').
            %
            % Parameters
            %   region_size : Size of the region (if shape is not 'global').
            %   position : Position offset (optional).
            if nargin < 2
                error('setPhase requires at least a value argument.');
            end
            shape = 'global';
            region_size = [];
            position = [];
            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);
            if ~strcmp(type, 'global')
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end
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
            % Adds a value to the phase of the field.
            %
            % Required
            %   value : Value to add.
            %
            % Optionals
            %   shape : Shape specification (default 'global').
            %
            % Parameters
            %   region_size : Size of the region (if shape is not 'global').
            %   position : Position offset (optional).
            if nargin < 2
                error('addPhase requires at least a value argument.');
            end
            shape = 'global';
            region_size = [];
            position = [];
            if nargin > 2, shape = varargin{1}; end
            type = obj.standardizeType(shape);
            if ~strcmp(type, 'global')
                if nargin < 4
                    error('Size is required for shape: %s', type);
                end
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
            % Maps alias strings to a canonical shape name.
            %
            % Required
            %   inputType : A string representing the shape.
            %
            % Returns
            %   type : Canonical shape name ('global', 'rectangle', 'circle', or 'annulus').
            typeMap = struct(...
                'global', 'global', 'g', 'global', 'glob', 'global', ...
                'rect', 'rectangle', 'rectangle', 'rectangle', 'r', 'rectangle', 'rectangular', 'rectangle', ...
                'circ', 'circle', 'circle', 'circle', 'c', 'circle', 'circular', 'circle', ...
                'annu', 'annulus', 'annulus', 'annulus', 'a', 'annulus', 'annular', 'annulus');
            if ~ischar(inputType) && ~isstring(inputType)
                error('Input type must be a string.');
            end
            inputType = lower(char(inputType));
            if isfield(typeMap, inputType)
                type = typeMap.(inputType);
            else
                error('Invalid type: %s. Must be ''global'', ''rectangle'', ''circle'', or ''annulus''.', inputType);
            end
        end
        
        function W = wdf(obj)
            % Displays the Wigner Distribution Function for the field.
            %
            % Required
            %   obj : Field object.
            field.disp(cross=true);
            [~, crossdata] = field.getPhaseCross();
            [d, f, t] = wvd(crossdata, 1000000, 'smoothedPseudo', MinThreshold=0, NumFrequencyPoints=1000);
            figure; imagesc(t, f, d(1:100, :));
            ylabel('Local Spatial Frequency (m^{-1})');
            xlabel('x (m)');
        end
        
        function obj = fft(obj)
            % Applies a 2D FFT to the field and updates the field.
            %
            % Required
            %   obj : Field object.
            obj = obj.setComplexField(fftshift(fft2(obj.complexField)));
        end
        
        function [xAxis, amplitudeData] = getAmplitudeCross(obj, axisType, pos)
            % Returns the amplitude cross section of a 2D field.
            %
            % Required
            %   obj : Field object.
            %
            % Optionals
            %   axisType : Axis for the cross section ('x' or 'y'; default: 'x').
            %   pos : Position index (default: middle of the field along the chosen axis).
            %
            % Returns
            %   xAxis : Spatial axis values.
            %   amplitudeData : Amplitude values along the cross section.
            if obj.dim ~= 2
                error('getAmplitudeCross is defined only for 2D fields.');
            end
            if nargin < 2 || isempty(axisType)
                axisType = 'x';
            end
            if nargin < 3 || isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(obj.size(1) / 2);
                else
                    pos = round(obj.size(2) / 2);
                end
            end
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
            % Returns the phase cross section of a 2D field.
            %
            % Required
            %   obj : Field object.
            %
            % Optionals
            %   axisType : Axis for the cross section ('x' or 'y'; default: 'x').
            %   pos : Position index (default: middle of the field along the chosen axis).
            %
            % Returns
            %   xAxis : Spatial axis values.
            %   phaseData : Phase values along the cross section.
            if obj.dim ~= 2
                error('getPhaseCross is defined only for 2D fields.');
            end
            if nargin < 2 || isempty(axisType)
                axisType = 'x';
            end
            if nargin < 3 || isempty(pos)
                if strcmp(axisType, 'x')
                    pos = round(obj.size(1) / 2);
                else
                    pos = round(obj.size(2) / 2);
                end
            end
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
            % Unwraps the phase of the field using a phase unwrapping utility.
            %
            % Required
            %   obj : Field object.
            obj.phase = osf.utils.phaseUnwrap(obj.phase);
        end
        
        function obj = show(obj, varargin)
            % Displays the field using a dashboard.
            %
            % Optionals
            %   'plotType' : Plot type ('default', 'amplitude', or 'phase'; default: 'default').
            %   'unwrap' : Logical flag to unwrap phase (default: true).
            %
            % Example
            %   obj.show('plotType', 'phase');
            p = inputParser;
            addOptional(p, 'plotType', 'default', @(x) ismember(x, {'default', 'amplitude', 'a', 'phase', 'p'}));
            addParameter(p, 'unwrap', true, @islogical);
            p.KeepUnmatched = true;
            parse(p, varargin{:});
            plotType = p.Results.plotType;
            unwrapFlag = p.Results.unwrap;
            
            phaseType = 'phase';
            if unwrapFlag
                phaseType = 'phase.unwrap';
            end
            
            if ismember(plotType, {'default'})
                osf.Dashboard([2,1], p.Unmatched).show(obj, 'amplitude').show(obj, phaseType);
            elseif ismember(plotType, {'amplitude', 'a'})
                osf.show(obj, 'amplitude', p.Unmatched);
            elseif ismember(plotType, {'phase', 'p'})
                osf.show(obj, 'phase', p.Unmatched);
            else
                error('Unsupported plotType %s', plotType);
            end
        end
        
        function obj = img(obj, path, varargin)
            % Imports an image and assigns it to the field.
            %
            % Required
            %   path : Path to the image file.
            %
            % Optionals
            %   'type' : Field type to update ('phase' or 'amplitude'; default: 'phase').
            %   'range' : Normalization range as a two-element vector (default: [0 1]).
            p = inputParser;
            addOptional(p, 'type', 'phase', @(x) ismember(x, {'phase', 'amplitude'}));
            addOptional(p, 'range', [0 1], @(x) isnumeric(x) && numel(x) == 2);
            parse(p, varargin{:});
            
            type = p.Results.type;
            range = p.Results.range;
            
            img = imread(path);
            if size(img, 3) == 3
                img = rgb2gray(img);
            end
            img = double(img);
            img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
            img = range(1) + img * (range(2) - range(1));
            targetSize = obj.size;
            img = imresize(img, targetSize, 'bilinear');
            
            if strcmp(type, 'phase')
                obj.phase = img;
            else
                obj.amplitude = img;
            end
        end
        
        function operationCheck(a, b)
            % Checks compatibility for operations between two Field objects.
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
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
            % Adds two Field objects element-wise.
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
            %
            % Returns
            %   result : Field object with summed amplitude and phase.
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude + b.amplitude;
            result.phase = a.phase + b.phase;
        end
        
        function result = minus(a, b)
            % Subtracts one Field object from another element-wise.
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
            %
            % Returns
            %   result : Field object with subtracted amplitude and phase.
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude - b.amplitude;
            result.phase = a.phase - b.phase;
        end
        
        function result = times(a, b)
            % Multiplies two Field objects element-wise.
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
            %
            % Returns
            %   result : Field object with multiplied amplitude and phase.
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude .* b.amplitude;
            result.phase = a.phase .* b.phase;
        end
        
        function result = mtimes(a, b)
            % Multiplies two Field objects element-wise (using * operator).
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
            %
            % Returns
            %   result : Field object with multiplied amplitude and phase.
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude .* b.amplitude;
            result.phase = a.phase .* b.phase;
        end
        
        function result = rdivide(a, b)
            % Divides one Field object by another element-wise.
            %
            % Required
            %   a : First Field object.
            %   b : Second Field object.
            %
            % Returns
            %   result : Field object with divided amplitude and phase.
            a.operationCheck(b);
            result = a;
            result.amplitude = a.amplitude ./ b.amplitude;
            result.phase = a.phase ./ b.phase;
        end
        
        function val = get.intensity(obj)
            % Returns the intensity of the field computed as the squared magnitude.
            %
            % Required
            %   obj : Field object.
            %
            % Returns
            %   val : Intensity.
            val = abs(obj.complexField) .^ 2;
        end
        
        function result = get.complexField(obj)
            % Returns the complex field computed as amplitude .* exp(1i*phase).
            %
            % Required
            %   obj : Field object.
            %
            % Returns
            %   result : Complex field.
            result = obj.getComplexField();
        end
    end
end
