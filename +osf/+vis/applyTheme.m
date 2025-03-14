function applyTheme(fig)
    if nargin < 1
        fig = gcf;
    end

    fig.Color = 'black';

    ax = findall(fig, 'Type', 'axes');
    for i = 1:length(ax)
        ax(i).XColor = 'white';
        ax(i).YColor = 'white';
        ax(i).ZColor = 'white';
        ax(i).Color = 'black';
        ax(i).GridColor = [0.5, 0.5, 0.5];
        ax(i).LineWidth = 1.2;
    end

    textObjs = findall(fig, 'Type', 'text');
    for i = 1:length(textObjs)
        textObjs(i).Color = 'white';
        textObjs(i).FontWeight = 'bold';
    end

    sgObjs = findall(fig, 'Tag', 'theme');
    for i = 1:length(sgObjs)
        sgObjs(i).Color = 'white';
    end

    colorbars = findall(fig, 'Type', 'colorbar');
    for i = 1:length(colorbars)
        colorbars(i).Color = 'white';
    end

    legends = findall(fig, 'Type', 'legend');
    for i = 1:length(legends)
        legends(i).TextColor = 'white';
        legends(i).Color = 'black';
    end
end
