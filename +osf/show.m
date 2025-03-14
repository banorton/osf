function fig = show(obj, varargin)

    if nargin < 1
        error('The first argument (obj) is required.');
    end

    objClass = class(obj);

    validClasses = {'osf.Sim', 'osf.Field'};
    if ~ismember(objClass, validClasses)
        error('Invalid object class: %s.', objClass);
    end

    switch objClass
    case 'osf.Sim'
        fig = osf.vis.sim(obj, varargin{:});
    case 'osf.Field'
        fig = osf.vis.field(obj, varargin{:});
    otherwise
        error('Support for object class %s has not been implemented.', objClass);
    end

end
