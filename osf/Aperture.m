classdef Aperture < Element
    properties
        name           % Name of the optical element
        phaseFunction  % Function handle for phase modification (if applicable)
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element ('1D' or '2D')
    end

    methods

        function obj = Aperture(varargin)
            obj.elementType = 'aperture';
            obj.apertureType = 'none';
            obj.apertureParams = struct();

            p = inputParser;
            addParameter(p, 'circ', NaN, @isnumeric);
            addParameter(p, 'rect', NaN, @isnumeric);
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            parse(p, varargin{:});

            obj.name = p.Results.name;
            obj.dim = p.Results.dim;

            if ~isnan(p.Results.circ) && ~isnan(p.Results.rect)
                error('Can not have ''circ'' and ''rect'' for aperture.');
            elseif ~isnan(p.Results.circ)
                obj = obj.addCircAperture(p.Results.circ);
            elseif ~isnan(p.Results.rect)
                obj = obj.addRectAperture(p.Results.rect);
            end

            if obj.dim == 2
                obj.phaseFunction = @(x, y, lambda) zeros(size(x));
            elseif obj.dim == 1
                obj.phaseFunction = @(x, lambda) zeros(size(x));
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function field = apply(obj, field)
            if ~strcmp(obj.apertureType, 'none')
                field = obj.applyAperture(field);
            end
        end

    end
end
