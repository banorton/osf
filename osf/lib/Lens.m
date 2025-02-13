classdef Lens < Element
    properties
        name           % Name of the optical element
        phaseFunction  % Function handle for phase modification (if applicable)
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)

        focalLength    % focalLength: Focal length in meters.
    end

    methods

        function obj = Lens(focalLength, varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'circ', 0, @isnumeric);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addRequired(p, 'focalLength', @isnumeric);
            parse(p, focalLength, varargin{:});

            obj.elementType = 'lens';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            if focalLength == 0
                error('Lens focal length must be non-zero.');
            end
            obj.focalLength = p.Results.focalLength;

            if obj.dim == 2
                obj.phaseFunction = @(x, y, lambda) (-2 * pi) / (2 * lambda * obj.focalLength) * (x.^2 + y.^2);
            elseif obj.dim == 1
                obj.phaseFunction = @(x, lambda) (-2 * pi) / (2 * lambda * obj.focalLength) * (x.^2);
            else
                error('Dimensionality must be either 1 or 2.');
            end

            if p.Results.circ > 0
                obj.addCircAperture(p.Results.circ);
            end

        end

        function field = apply(obj, field)
            distances = obj.getDistanceMatrix(field);

            if ~isempty(obj.phaseFunction)
                if obj.dim == 2
                    phaseShift = obj.phaseFunction(distances.X, distances.Y, field.lambda);
                elseif obj.dim == 1
                    phaseShift = obj.phaseFunction(distances.X, field.lambda);
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
