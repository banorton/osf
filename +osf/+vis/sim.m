function fig = sim(sim, varargin)

    p = inputParser;
    addRequired(p, 'sim');
    addOptional(p, 'plotType', 'simDefault', @(x) ischar(x) && ismember(x, {'default', 'simDefault'}));
    addParameter(p, 'title', '', @(x) ischar(x) || isstring(x));

    p.KeepUnmatched = true;
    parse(p, sim, varargin{:});

    plotType = p.Results.plotType;
    titleStr = p.Results.title;

    if strcmp(plotType, 'default')
        plotType = 'simDefault';
    end

    fig = figure('Position', [275 400 1300 350], 'Color', 'white');
    ax = subplot(1,1,1);

    if ~isempty(titleStr)
        sgtitle(titleStr, 'FontWeight', 'bold', 'FontSize', 18, 'Tag', 'theme');
    end

    data = osf.vis.prepData(sim, plotType, p.Unmatched);
    osf.vis.attachData(ax, data);
    osf.vis.applyTheme(fig);

end
