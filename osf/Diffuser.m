classdef Diffuser < Element
    properties
        name           % Name of the optical element
        phaseFunction  % Function handle for phase modification (if applicable)
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)

        roughness           % RMS roughness in meters
        correlation_length  % Correlation length in meters
    end

    methods

        function obj = Diffuser(roughness, correlation_length, varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addRequired(p, 'roughness', @isnumeric);
            addRequired(p, 'correlation_length', @isnumeric);
            parse(p, roughness, correlation_length, varargin{:});

            obj.elementType = 'diffuser';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            obj.roughness = roughness;
            obj.correlation_length = correlation_length;
            obj.phaseFunction = @(sz, res, lambda) imgaussfilt(2 * pi / lambda * obj.roughness * randn(sz), obj.correlation_length / res);
        end

        function field = apply(obj, field)
            if ~isempty(obj.phaseFunction)
                if obj.dim == 2
                    phaseShift = obj.phaseFunction(size(field.phase), field.resolution, field.lambda);
                elseif obj.dim == 1
                    phaseShift = obj.phaseFunction(size(field.phase), field.resolution, field.lambda);
                else
                    error('Unsupported dimensionality in element.');
                end
                field.phase = obj.wrap(field.phase + phaseShift);
            end

            if ~strcmp(obj.apertureType, 'none')
                field = obj.applyAperture(field);
            end
        end

    end
end
