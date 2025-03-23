function attachData(ax, data)
    cla(ax);

    switch data.meta.plotType
    case {'osf.Field.default', 'a', 'amplitude', 'phase', 'p', 'phase.unwrap', 'p.unwrap'}
        imagesc(data.xAxis, data.yAxis, data.imageData);
        colormap(gca, data.cmap);
        colorbar;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        axis equal; axis tight;

    case {'cross', 'amplitude.cross', 'a.cross', 'phase.cross', 'p.cross'}
        hold on;
        plot(data.xAxis, data.xCross, 'LineWidth', 1.5, 'DisplayName', 'X Cross Section');
        plot(data.yAxis, data.yCross, 'LineWidth', 1.5, 'DisplayName', 'Y Cross Section');
        hold off;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        grid off;
        legend('show');

    case {'wdf', 'amplitude.wdf', 'a.wdf', 'phase.wdf', 'p.wdf'}
        imagesc(data.xAxis, data.yAxis, data.imageData);
        axis xy;
        colormap(gca, data.colormap);
        colorbar;
        title(data.colorbarLabel);
        xlabel(data.xlabel);
        ylabel(data.ylabel);

    case {'osf.Sim.default'}
        attachSim(ax, data);

    case {'double.default', 'single.default', 'uint8.default', 'uint16.default'}
        imagesc(data.xAxis, data.yAxis, data.imageData);
        colormap(gca, data.cmap);
        colorbar;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        axis image;

    otherwise
        error('Unhandled plotType: %s', data.meta.plotType);
    end

    pause(.0001);

    function attachSim(ax, data)
        axes(ax);
        hold(ax, 'on');

        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';

        for i = 1:length(data.elements)
            element = data.elements(i);

            switch element.type
            case 'lens'
                plotLens(ax, element.position, data.elementHeight);
            case 'diffuser'
                plotDiffuser(ax, element.position, data.elementHeight);
            case 'plane'
                plotPlane(ax, element.position, data.elementHeight);
            otherwise
                warning('Unknown element type: %s', element.type);
                plotUnknown(ax, element.position, data.elementHeight);
            end

            text(ax, element.position - 0.000, data.elementHeight / 2 + data.labelOffset, ...
            element.name, 'HorizontalAlignment', 'center', 'FontSize', 12, ...
            'FontWeight', 'bold', 'Rotation', 45);
        end

        xlim(ax, data.axisLimits(1:2));
        ylim(ax, data.axisLimits(3:4));

        plot(ax, [min(data.componentCenters), max(data.componentCenters)], ...
        [data.baselineY, data.baselineY], 'k', 'LineWidth', 2, 'Tag', 'theme');

        for i = 1:length(data.componentCenters)
            xTick = data.componentCenters(i);
            plot(ax, [xTick, xTick], [data.baselineY, data.baselineY + data.tickHeight], 'k', 'LineWidth', 2, 'Tag', 'theme');

            if i < length(data.componentCenters)
                midX = (data.componentCenters(i) + data.componentCenters(i+1)) / 2;
                text(ax, midX, data.baselineY - 0.01, sprintf('%.0fmm', data.distances(i+1)*1000), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            end
        end

        % --- Paraxial Ray Overlay ---
        plot(ax, data.paraxial.dist, data.paraxial.marginalRay, 'r', 'LineWidth', 0.5);
        plot(ax, data.paraxial.dist, -data.paraxial.marginalRay, 'r', 'LineWidth', 0.5);
        plot(ax, data.paraxial.dist, data.paraxial.chiefRay, 'b', 'LineWidth', 0.5);
        plot(ax, data.paraxial.dist, -data.paraxial.chiefRay, 'b', 'LineWidth', 0.5);

        hold(ax, 'off');

        function plotLens(ax, x, elementHeight)
            width = 0.005;
            t = linspace(0, pi, 20);
            X = [cos(t), -cos(t)] * width/2 + x;
            Y = [sin(t), -sin(t)] * elementHeight/2;
            fill(ax, X, Y, [.8 .8 .8], 'EdgeColor', 'k', 'LineWidth', 2);
        end

        function plotDiffuser(ax, x, elementHeight)
            width = 0.002;
            rectangle(ax, 'Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
            'FaceColor', 'g', 'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        function plotPlane(ax, x, elementHeight)
            plot(ax, [x, x], [-elementHeight/2, elementHeight/2], 'k', 'LineWidth', 2, 'Tag', 'theme');
        end

        function plotUnknown(ax, x, elementHeight)
            width = 0.003;
            rectangle(ax, 'Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
            'FaceColor', 'r', 'EdgeColor', 'k', 'LineWidth', 1.5);
        end
    end

end
