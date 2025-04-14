classdef Diffuser < osf.elements.Element
    properties
        name           % Name of the optical element
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)
        id

        roughness          % RMS roughness in meters
        correlationLength  % Correlation length in meters
    end

    methods

        function obj = Diffuser(roughness, correlationLength, varargin)
            p = inputParser;
            addRequired(p, 'roughness', @isnumeric);
            addRequired(p, 'correlationLength', @isnumeric);
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'id', 0, isnumeric(x));
            parse(p, roughness, correlationLength, varargin{:});

            obj.name = obj.genName(p.Results.name);
            obj.dim = p.Results.dim;
            obj.id = p.Results.id;

            obj.elementType = 'diffuser';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.roughness = roughness;
            obj.correlationLength = correlationLength;
            obj.correlationLength = correlationLength;
        end

        function phaseShift = phaseFunction(obj, sz, res, lambda)
            if obj.roughness == 0 || obj.correlationLength == 0
                phaseShift = 0;
            else
                phaseShift = imgaussfilt(2*pi/lambda * obj.roughness * randn(sz), obj.correlationLength/res);
            end
        end

        function field = apply(obj, field)
            if obj.dim == 2
                phaseShift = obj.phaseFunction(size(field.phase), field.resolution, field.lambda);
            elseif obj.dim == 1
                phaseShift = obj.phaseFunction(size(field.phase), field.resolution, field.lambda);
            else
                error('Unsupported dimensionality in element.');
            end
            field.phase = field.phase + phaseShift;

            if ~strcmp(obj.apertureType, 'none')
                field = obj.applyAperture(field);
            end
        end

        function print(obj)
            % PRINT Prints simulation parameters.
            roughness_um           = obj.roughness * 1e6;
            correlation_length_um  = obj.correlationLength * 1e6;

            fprintf('\nDiffuser Parameters:\n');
            fprintf('-------------------------------\n');
            fprintf('  Surface Roughness:       %.2f um\n', roughness_um);
            fprintf('  Correlation Length:      %.2f um\n', correlation_length_um);
        end

    end
end
