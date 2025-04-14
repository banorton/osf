classdef Aperture < osf.elements.Element
    properties
        name           % Name of the optical element
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element ('1D' or '2D')
        id
    end

    methods

        function obj = Aperture(varargin)
            p = inputParser;
            addParameter(p, 'circ', NaN, @isnumeric);
            addParameter(p, 'rect', NaN, @isnumeric);
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'id', 0, @isnumeric);
            parse(p, varargin{:});

            obj.dim = p.Results.dim;
            obj.id = p.Results.id;
            obj.elementType = 'aperture';
            obj.apertureType = 'none';
            obj.apertureParams = struct();

            obj.name = obj.genName(p.Results.name);

            if ~isnan(p.Results.circ) && ~isnan(p.Results.rect)
                error('Can not have ''circ'' and ''rect'' for aperture.');
            elseif ~isnan(p.Results.circ)
                obj = obj.addCircAperture(p.Results.circ);
            elseif ~isnan(p.Results.rect)
                obj = obj.addRectAperture(p.Results.rect);
            end

        end

        function phaseShift = phaseFunction(obj)
            phaseShift = 0;
        end

        function field = apply(obj, field)
            if ~strcmp(obj.apertureType, 'none')
                field = obj.applyAperture(field);
            end
        end

    end
end
