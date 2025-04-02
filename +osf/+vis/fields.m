function fig = fields(fArray, varargin)
    if iscell(fArray)
        if ~all(cellfun(@(f) isa(f, 'osf.Field'), fArray))
            error('All elements of the cell array must be osf.Field objects.');
        end
        fieldsFlat = [fArray{:}];
    elseif isa(fArray, 'osf.Field')
        fieldsFlat = fArray;
    else
        error('Input must be an array or cell array of osf.Field objects.');
    end

    n = numel(fieldsFlat);
    fig = figure('Position', [300 300 300*n 600]);

    for i = 1:n
        % Amplitude subplot
        ax1 = subplot(2, n, i);
        dataA = osf.vis.prepData(fieldsFlat(i), 'amplitude', varargin{:});
        osf.vis.attachData(ax1, dataA);

        % Phase subplot
        ax2 = subplot(2, n, n+i);
        dataP = osf.vis.prepData(fieldsFlat(i), 'phase', 'unwrap', true, varargin{:});
        osf.vis.attachData(ax2, dataP);
    end

    osf.vis.applyTheme(fig);
end
