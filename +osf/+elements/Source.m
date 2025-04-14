classdef Source < osf.elements.Element
    properties
        name
        elementType
        apertureType
        apertureParams
        dim
        id

        wavelength
    end

    methods
        function obj = Source(wavelength, varargin)
            p = inputParser;
            addOptional(p, 'wavelength', 632e-9, @isnumeric);
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'id', 0, @isnumeric);
            parse(p, wavelength, varargin{:});

            obj.dim = p.Results.dim;
            obj.wavelength = p.Results.wavelength;
            obj.id = p.Results.id;

            obj.elementType = 'source';
            obj.apertureType = 'none';
            obj.apertureParams = struct();

            obj.name = obj.genName(p.Results.name);
        end

        function field = apply(obj, field)
            field = field.newField().setPhase(1) * field;
        end

        function field = phaseFunction(obj, varargin)
        end
    end
end
