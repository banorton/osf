classdef Dashboard < handle
    properties
        subplotRows
        subplotCols
        numSubplots
        occupiedSubplots
        fig
        title
        figPosition
    end

    methods
        function obj = Dashboard(subplotSize, varargin)
            if numel(subplotSize) ~= 2
                error('subplotSize must be a two-element vector: [rows, cols]');
            end

            p = inputParser;
            addParameter(p, 'figPosition', [750 350 450 450], @(x) isnumeric(x) && numel(x) == 4);
            addParameter(p, 'title', '', @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});

            obj.figPosition = p.Results.figPosition;
            obj.title = p.Results.title;

            obj.subplotRows = subplotSize(1);
            obj.subplotCols = subplotSize(2);
            obj.numSubplots = obj.subplotRows * obj.subplotCols;
            obj.occupiedSubplots = false(1, obj.numSubplots);
            obj.fig = figure('Position', obj.figPosition);

            if ~isempty(obj.title)
                sgtitle(obj.fig, obj.title, 'FontWeight', 'bold', 'FontSize', 20, 'Tag', 'theme');
            end
        end

        function obj = show(obj, varargin)
            if nargin == 1
                obj.reset();
                return;
            end

            if nargin < 3 || isempty(varargin{2})
                plotType = 'default';
                dataObj = varargin{1};
                extraArgs = {};
            else
                dataObj = varargin{1};
                plotType = varargin{2};
                extraArgs = varargin(3:end);
            end

            nextIndex = find(~obj.occupiedSubplots, 1);
            if isempty(nextIndex)
                error('All subplots are filled. Cannot add more plots.');
            end

            obj.occupiedSubplots(nextIndex) = true;
            subplot(obj.subplotRows, obj.subplotCols, nextIndex);

            data = osf.vis.prepData(dataObj, plotType, extraArgs{:});
            osf.vis.attachData(gca, data);
            osf.vis.applyTheme(obj.fig);
        end

        function obj = showOn(obj, subplotPos, dataObj, plotType, varargin)
            if nargin < 3
                error('Not enough input arguments. Expected at least (subplotPos, dataObj).');
            end

            if isnumeric(subplotPos) && isscalar(subplotPos)
                subplotIndex = subplotPos;
            elseif isnumeric(subplotPos) && numel(subplotPos) == 2
                row = subplotPos(1);
                col = subplotPos(2);
                subplotIndex = (row - 1) * obj.subplotCols + col;
            else
                error('Invalid subplot position. Provide a scalar index or a [row, col] pair.');
            end

            if subplotIndex > obj.numSubplots || subplotIndex < 1
                error('Invalid subplot index: %d. Must be between 1 and %d.', subplotIndex, obj.numSubplots);
            end

            subplot(obj.subplotRows, obj.subplotCols, subplotIndex);

            if nargin < 4 || isempty(plotType)
                plotType = 'default';
            end

            data = osf.vis.prepData(dataObj, plotType, varargin{:});
            osf.vis.attachData(gca, data);
            osf.vis.applyTheme(obj.fig);
        end

        function reset(obj)
            if ~isvalid(obj.fig) || ~strcmp(obj.fig.Type, 'figure')
                obj.fig = figure('Position', obj.figPosition);
                if ~isempty(obj.title)
                    sgtitle(obj.fig, obj.title, 'FontWeight', 'bold', 'FontSize', 20, 'Tag', 'theme');
                end
            else
                figure(obj.fig);
            end
        end
    end
end
