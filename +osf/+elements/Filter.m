classdef Filter < osf.elements.Element
    properties
        name           % Name of the optical element
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', '1D-aperture')
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)
        id

        field          % The field transformation filter applies
        operation      % Operation mode ('mult' or 'add')
    end

    methods
        function obj = Filter(field, varargin)
            % Constructor for Filter class
            p = inputParser;
            addRequired(p, 'field', @(x) isa(x, 'osf.Field'));
            addParameter(p, 'name', 'Filter', @(x) ischar(x) || isstring(x));
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'operation', 'mult', @(x) ischar(x) || isstring(x) && ismember(x, {'mult', 'add'}));
            addParameter(p, 'id', 0, @isnumeric(x));

            parse(p, field, varargin{:});

            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            obj.operation = char(p.Results.operation);
            obj.id = p.Results.id;
            obj.field = field;
            obj.elementType = 'filter';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
        end

        function phaseShift = phaseFunction(~)
            % No specific phase function for a general filter
            phaseShift = 0;
        end

        function field = apply(obj, field)
            switch obj.operation
                case 'mult'
                    field.amplitude = field.amplitude .* obj.field.amplitude;
                    field.phase = obj.field.phase .* field.phase;

                case 'add'
                    field.amplitude = field.amplitude + obj.field.amplitude;
                    field.phase = field.phase + obj.field.phase;
            end
        end
    end
end
