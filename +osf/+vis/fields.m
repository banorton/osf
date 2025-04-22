function fig = fields(fArray, varargin)
    fields = [fArray{:}];

    n = numel(fields);
    fig = figure('Position', [300 300 300*n 600]);

    for i = 1:n
        ax1 = subplot(2, n, i);
        dataA = osf.vis.prepData(fields(i), 'amplitude', varargin{:});
        osf.vis.attachData(ax1, dataA);

        ax2 = subplot(2, n, n+i);
        dataP = osf.vis.prepData(fields(i), 'phase', varargin{:});
        osf.vis.attachData(ax2, dataP);
    end

    osf.vis.applyTheme(fig);
end
