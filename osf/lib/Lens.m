classdef Lens < Element
    properties
        name           % Name of the optical element
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

            if p.Results.circ > 0
                obj.addCircAperture(p.Results.circ);
            end

        end

        function phaseShift = phaseFunction(obj, varargin)
            % PHASEFUNCTION Computes the phase shift for the element.
            %   For 2D fields, call as:
            %       phaseShift = obj.phaseFunction(x, y, lambda)
            %   For 1D fields, call as:
            %       phaseShift = obj.phaseFunction(x, lambda)

            if obj.dim == 2
                if numel(varargin) < 3
                    error('For 2D, provide x, y, and lambda.');
                end
                x = varargin{1};
                y = varargin{2};
                lambda = varargin{3};
                phaseShift = (-2 * pi) / (2 * lambda * obj.focalLength) * (x.^2 + y.^2);
            elseif obj.dim == 1
                if numel(varargin) < 2
                    error('For 1D, provide x and lambda.');
                end
                x = varargin{1};
                lambda = varargin{2};
                phaseShift = (-2 * pi) / (2 * lambda * obj.focalLength) * (x.^2);
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function field = apply(obj, field)
            distances = obj.getDistanceMatrix(field);

            if obj.dim == 2
                phaseShift = obj.phaseFunction(distances.X, distances.Y, field.lambda);
            elseif obj.dim == 1
                phaseShift = obj.phaseFunction(distances.X, field.lambda);
            else
                error('Unsupported dimensionality in element.');
            end
            field.phase = obj.wrap(field.phase + phaseShift);

            if ~strcmp(obj.apertureType, 'none')
                field = obj.applyAperture(field);
            end
        end

    end
end
