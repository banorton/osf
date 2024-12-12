classdef Plane < Element
    properties
    end

    methods

        function obj = Plane(varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            parse(p, varargin{:});

            obj.name = p.Results.name;
            obj.dim = p.Results.dim;

            obj.phaseFunction = @() 0;
        end

        function field = apply(obj, field)
        end

    end
end
