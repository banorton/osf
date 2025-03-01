classdef Filter < osf.Element
    properties
        name           % Name of the optical element
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)

        type           % Type of filter ('ampltiude', or 'phase')
        shape          % Shape of filter ('rect', 'circ', etc.)
    end

    methods

        function obj = Filter(varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            parse(p, varargin{:});

            obj.elementType = 'filter';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
        end

        function phaseShift = phaseFunction(obj)
            phaseShift = 0;
        end

        function field = apply(obj, field)

        end

    end
end
