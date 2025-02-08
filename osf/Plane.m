classdef Plane < Element
    properties
        name           % Name of the optical element
        phaseFunction  % Function handle for phase modification (if applicable)
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)
    end

    methods

        function obj = Plane(varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            parse(p, varargin{:});

            obj.elementType = 'plane';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;

            if obj.dim == 2
                obj.phaseFunction = @(x, y, lambda) zeros(size(x));
            elseif obj.dim == 1
                obj.phaseFunction = @(x, lambda) zeros(size(x));
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function field = apply(obj, field)
        end

    end
end
