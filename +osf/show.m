function fig = show(obj, varargin)

    if nargin < 1
        error('The first argument (obj) is required.');
    end

    objClass = class(obj);

    if isa(obj, 'osf.Sim')
        fig = osf.vis.sim(obj, varargin{:});
    elseif isa(obj, 'osf.Field')
        fig = osf.vis.field(obj, varargin{:});
    elseif isnumeric(obj) && isreal(obj)
        fig = osf.vis.matrix(obj, varargin{:});
    else
        error('Unsupported input type: %s.', objClass);
    end

end
