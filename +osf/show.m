function fig = show(obj, varargin)
% Displays an object using the appropriate visualization function.
%
% Required
%   obj : The object to display. It can be one of the following:
%         - An osf.Sim object
%         - An osf.Field object (scalar or array)
%         - A cell array of osf.Field objects
%         - A real numeric matrix
%
% Optionals
%   varargin : Additional parameters passed to the visualization function.
%
% Example
%   fig = show(simObj, 'optionName', optionValue);

    if nargin < 1
        error('The first argument (obj) is required.');
    end

    if isa(obj, 'osf.Sim')
        fig = osf.vis.sim(obj, varargin{:});
    elseif isa(obj, 'osf.Field')
        if isscalar(obj)
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
