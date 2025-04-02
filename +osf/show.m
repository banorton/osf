function fig = show(obj, varargin)
    if nargin < 1
        error('The first argument (obj) is required.');
    end

    if isa(obj, 'osf.Sim')
        fig = osf.vis.sim(obj, varargin{:});
    elseif isa(obj, 'osf.Field')
        if numel(obj) == 1
            fig = osf.vis.field(obj, varargin{:});
        else
            fig = osf.vis.fields(obj, varargin{:});
        end
    elseif iscell(obj) && all(cellfun(@(f) isa(f, 'osf.Field'), obj))
        fig = osf.vis.fields(obj, varargin{:});
    elseif isnumeric(obj) && isreal(obj)
        fig = osf.vis.matrix(obj, varargin{:});
    else
        error('Unsupported input type: %s.', class(obj));
    end
end
