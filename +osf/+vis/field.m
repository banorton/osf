function fig = field(f, varargin)

    p = inputParser;
    addRequired(p, 'f');
    addOptional(p, 'plotType', 'default', @(x) ischar(x) && ismember(x, ...
        {'default', 'a', 'amplitude', 'p', 'phase', 'cross', 'amplitude.cross', 'a.cross', ...
        'phase.cross', 'p.cross', 'wdf', 'amplitude.wdf', 'a.wdf', ...
        'phase.wdf', 'p.wdf', 'fft'}));
    addParameter(p, 'title', '', @(x) ischar(x) || isstring(x));

    p.KeepUnmatched = true;
    parse(p, f, varargin{:});
    plotType = p.Results.plotType;
    titleStr = p.Results.title;

    fig = figure('Position', [750 350 450 450]);
    ax = subplot(1,1,1);

    if ~isempty(titleStr)
        sgtitle(titleStr, 'FontWeight', 'bold', 'FontSize', 18, 'Tag', 'theme');
    end

    data = osf.vis.prepData(f, plotType, p.Unmatched);
    osf.vis.attachData(ax, data);
    osf.vis.applyTheme(fig);

end
