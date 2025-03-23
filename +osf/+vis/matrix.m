function fig = matrix(m, varargin)

    p = inputParser;
    addRequired(p, 'm');
    addParameter(p, 'title', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'colormap', 'gray', @(x) ischar(x) || isstring(x));
    p.KeepUnmatched = true;
    parse(p, m, varargin{:});

    titleStr = p.Results.title;
    cmap = p.Results.colormap;

    fig = figure('Position', [750 350 450 450]);
    ax = subplot(1,1,1);

    if ~isempty(titleStr)
        sgtitle(titleStr, 'FontWeight', 'bold', 'FontSize', 18, 'Tag', 'theme');
    end

    data = osf.vis.prepData(m, sprintf("%s.%s", class(m), "default"), p.Unmatched);
    osf.vis.attachData(ax, data);
    osf.vis.applyTheme(fig);

end
