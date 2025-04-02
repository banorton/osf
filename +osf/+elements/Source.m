classdef Source < osf.elements.Element
    properties
        name
        elementType
        apertureType
        apertureParams
        dim

        wavelength
    end

    methods
        function obj = Source(wavelength, varargin)
            p = inputParser;
            addOptional(p, 'wavelength', 632e-9, @isnumeric);
            addParameter(p, 'name', 'Source', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            parse(p, wavelength, varargin{:});

            obj.elementType = 'source';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            obj.wavelength = p.Results.wavelength;
        end

        function field = apply(obj, field)
            field.addAmplitude(1);
        end

        function field = phaseFunction(obj, varargin)
        end
    end
end
