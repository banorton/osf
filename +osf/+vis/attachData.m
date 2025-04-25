function attachData(ax, data)
    axes(ax);

    switch data.meta.plotType
    case {'osf.Field.default', 'a', 'amplitude', 'phase', 'p'}
        data.imageData = flipud(data.imageData);
        imagesc(data.xAxis, data.yAxis, data.imageData);
        colormap(ax, data.cmap);
        colorbar;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        axis xy;

    case {'cross', 'amplitude.cross', 'a.cross', 'phase.cross', 'p.cross'}
        hold on;
        if ismember(data.crossAxis, {'both', 'y'})
            plot(data.yAxis, data.yCross, '-', 'Color', [1 0 0 0.8], 'MarkerSize', 3, 'DisplayName', 'Y Cross Section');
        end
        if ismember(data.crossAxis, {'both', 'x'})
            plot(data.xAxis, data.xCross, '-', 'Color', [0.3 0.5 1 0.8], 'MarkerSize', 3, 'DisplayName', 'X Cross Section');
        end
        hold off;

        xCenter = mean(data.xAxis);
        xRange = (max(data.xAxis) - min(data.xAxis)) / data.zoom / 2;
        xlim([xCenter - xRange, xCenter + xRange]);

        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        box on;
        grid off;
        legend('show');

    case {'wdf', 'amplitude.wdf', 'a.wdf', 'phase.wdf', 'p.wdf'}
        imagesc(data.xAxis, data.yAxis, data.imageData);
        if  isfield(data, 'xBox')
            hold on;
            fill(data.xBox, data.yBox, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r');
        end
        axis xy;
        colormap(gca, data.colormap);
        c = colorbar;
        c.Label.String = data.colorbarLabel;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);

    case {'osf.Sim.default'}
        cla(ax);
        attachSim(ax, data);

    case {'double.default', 'single.default', 'uint8.default', 'uint16.default'}
        imagesc(data.xAxis, data.yAxis, data.imageData);
        colormap(gca, data.cmap);
        colorbar;
        title(data.title);
        xlabel(data.xlabel);
        ylabel(data.ylabel);
        % axis xy;

    otherwise
        error('Unhandled plotType: %s', data.meta.plotType);
    end

    pause(.0001);

    function attachSim(ax, data)
        axes(ax);
        hold(ax, 'on');

        % Hide axis lines
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';

        if data.isempty == true
            return
        end

        for i = 1:length(data.elements)
            element = data.elements(i);

            switch element.type
            case 'lens'
                plotLens(ax, element.position, data.elementHeight);
            case 'diffuser'
                plotDiffuser(ax, element.position, data.elementHeight);
            case 'plane'
                plotPlane(ax, element.position, data.elementHeight);
            case 'source'
                plotSource(ax, element.position, data.elementHeight);
            case 'detector'
                plotDetector(ax, element.position, data.elementHeight);
            case 'filter'
                plotFilter(ax, element.position, data.elementHeight);
            case 'aperture'
                plotAperture(ax, element.position, data.elementHeight);
            case 'grating'
                plotGrating(ax, element.position, data.elementHeight);
            otherwise
                plotUnknown(ax, element.position, data.elementHeight);
            end

            text(ax, element.position, data.elementHeight / 2 + data.labelOffset, ...
            element.name, 'HorizontalAlignment', 'center', 'FontSize', 12, ...
            'FontWeight', 'bold', 'Rotation', 45);
        end

        xlim(ax, data.axisLimits(1:2));
        ylim(ax, data.axisLimits(3:4));

        if length(data.elements) > 1
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
        end

        % Only overlay paraxial rays if more than one element exists.
        if length(data.elements) > 1 && isfield(data, 'paraxial') && ~isempty(data.paraxial)
            plot(ax, data.paraxial.dist, data.paraxial.marginalRay, 'r', 'LineWidth', 1.5);
            plot(ax, data.paraxial.dist, -data.paraxial.marginalRay, 'r', 'LineWidth', 1.5);
            plot(ax, data.paraxial.dist, data.paraxial.chiefRay, 'b', 'LineWidth', 1.5);
            plot(ax, data.paraxial.dist, -data.paraxial.chiefRay, 'b', 'LineWidth', 1.5);
        end

        hold(ax, 'off');

        function plotSource(ax, x, elementHeight)
            width = 0.002;
            rectangle(ax, 'Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
            'FaceColor', [0.2 0.6 1.0], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        function plotDetector(ax, x, elementHeight)
            width = 0.002;
            rectangle(ax, 'Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
            'FaceColor', [0.0 0.7 0.2], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        function plotFilter(ax, x, elementHeight)
            width = 0.002;
            rectangle(ax, 'Position', [x - width/2, -elementHeight/2, width, elementHeight], ...
            'FaceColor', [0.6 0.0 0.6], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        function plotLens(ax, x, elementHeight)
            width = 0.002;
            t = linspace(0, pi, 20);
            X = [cos(t), -cos(t)] * width/2 + x;
            Y = [sin(t), -sin(t)] * elementHeight/2;
            fill(ax, X, Y, [.8 .8 .8], 'EdgeColor', 'k', 'LineWidth', 2);
        end

        function plotGrating(ax, x, elementHeight)
            width = 0.002;
            numLines = 10;
            left = x - width/2;
            stripeWidth = width / numLines;
            yBottom = -elementHeight/2;
            for i = 0:numLines-1
                stripeX = left + i * stripeWidth;
                color = mod(i, 2) == 0;
                rectangle(ax, 'Position', [stripeX, yBottom, stripeWidth, elementHeight], ...
                'FaceColor', color * [1 1 1], 'EdgeColor', 'none');
            end
            rectangle(ax, 'Position', [left, yBottom, width, elementHeight], ...
            'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        function plotAperture(ax, x, elementHeight)
            width = 0.002;
            yBottom = -elementHeight/2;
            rectangle(ax, 'Position', [x - width/2, yBottom, width, elementHeight], ...
            'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1.5);
            radius = elementHeight * 0.2;
            theta = linspace(0, 2*pi, 100);
            fill(ax, x + radius*cos(theta), radius*sin(theta), 'w', 'EdgeColor', 'none');
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
