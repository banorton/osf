function applyTheme(fig)
    if nargin < 1
        fig = gcf;
    end

    fig.Color = 'black';

    ax = findall(fig, 'Type', 'axes');
    for i = 1:length(ax)
        if ~strcmp(ax(i).XColor, 'none')
            ax(i).XColor = 'white';
        end
        if ~strcmp(ax(i).YColor, 'none')
            ax(i).YColor = 'white';
        end
        if ~strcmp(ax(i).ZColor, 'none')
            ax(i).ZColor = 'white';
        end
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
