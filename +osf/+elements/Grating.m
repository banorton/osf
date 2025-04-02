classdef Grating < osf.elements.Element
    properties
        name           % Name of the optical element
        elementType    % Type of optical element (e.g., 'grating')
        apertureType   % Type of aperture ('none', 'circ', 'rect', etc.)
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)

        linesPerMM  % Lines per mm (scalar for 1D, [linesX, linesY] for 2D)
        gratingType        % Grating type ('sinusoidal' or 'blazed')
        blazeAngle  % Blaze angle (only for blazed grating)
        type        % Target field property to apply grating ('amplitude' or 'phase')
    end

    methods
        function obj = Grating(linesPerMM, varargin)
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addOptional(p, 'gratingType', 'sinusoidal', @(x) ismember(x, {'sinusoidal', 's', 'sine', 'blazed'}));
            addOptional(p, 'blazeAngle', 0, @(x) isnumeric(x) && x >= 0);
            addOptional(p, 'type', 'amplitude', @(x) ismember(x, {'amplitude', 'phase'}));
            addRequired(p, 'linesPerMM', @(x) isnumeric(x) && (isscalar(x) || numel(x) == 2));
            parse(p, linesPerMM, varargin{:});

            obj.dim = p.Results.dim;
            obj.name = 'Grating';
            obj.elementType = 'grating';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.type = p.Results.type;
            obj.linesPerMM = p.Results.linesPerMM;

            if ~ismember(obj.dim, [1, 2])
                error('Dimensionality must be either 1 or 2.');
            end

            obj.gratingType = obj.standardizeType(p.Results.gratingType);
            obj.blazeAngle = p.Results.blazeAngle;
        end

        function type = standardizeType(~, inputType)
            % Maps multiple aliases to a single canonical type
            typeMap = struct('sinusoidal', 'sinusoidal', 's', 'sinusoidal', 'sine', 'sinusoidal', 'blazed', 'blazed');
            if isfield(typeMap, inputType)
                type = typeMap.(inputType);
            else
                error('Invalid type: %s. Must be ''sinusoidal'' or ''blazed''.', inputType);
            end
        end

        function grating = createGrating(obj, fieldLength, resolution)
            % Generate grating field with specified size and resolution
            x = linspace(-fieldLength/2, fieldLength/2, round(fieldLength / resolution));
            if obj.dim == 1
                grating = obj.createGrating1D(x);
            else
                y = linspace(-fieldLength/2, fieldLength/2, round(fieldLength / resolution));
                [X, Y] = meshgrid(x, y);
                grating = obj.createGrating2D(X, Y);
            end
        end

        function grating = createGrating1D(obj, x)
            % Create 1D grating
            frequency = obj.linesPerMM * 1e3; % Convert to lines per meter
            if strcmp(obj.gratingType, 'sinusoidal')
                grating = 0.5 + 0.5 * sin(2 * pi * frequency * x);
            else % Blazed grating
                grating = mod(frequency * x, 1);
            end
        end

        function grating = createGrating2D(obj, X, Y)
            % Create 2D grating
            frequencyX = obj.linesPerMM(1) * 1e3; % Convert to lines per meter
            frequencyY = obj.linesPerMM(2) * 1e3;
            if strcmp(obj.gratingType, 'sinusoidal')
                grating = 0.5 + 0.5 * sin(2 * pi * (frequencyX * X + frequencyY * Y));
            else % Blazed grating
                grating = mod(frequencyX * X + frequencyY * Y, 1);
            end
        end

        function field = apply(obj, field)
            % Create the grating phase pattern
            gratingPattern = obj.createGrating(field.fieldLength, field.resolution);

            if strcmp(obj.type, 'phase')
                field.phase = field.phase + gratingPattern;
            else
                field.amplitude = field.amplitude .* gratingPattern;
            end
        end

        function phaseShift = phaseFunction(~)
            % No specific phase function for a general filter
            phaseShift = 0;
        end
    end
end
