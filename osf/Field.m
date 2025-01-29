classdef Field
    properties
        amplitude    % Amplitude of the field (1D or 2D array)
        phase        % Phase of the field (1D or 2D array)
        fieldLength  % Length of the field in meters (scalar for 1D, [width, height] for 2D)
        resolution   % Spatial resolution in meters per pixel (scalar)
        lambda       % Wavelength of light in meters
        dim          % Dimensionality of the field (1 or 2)
    end

    methods
        function obj = Field(dim, fieldLength, resolution, lambda)
            % Validate dimensionality
            if ~ismember(dim, [1, 2])
                error('Dimensionality must be either 1 or 2.');
            end

            obj.dim = dim;
            obj.fieldLength = fieldLength;
            obj.resolution = resolution;
            obj.lambda = lambda;

            % Initialize amplitude and phase based on dimensionality
            if dim == 1
                samples = round(fieldLength / resolution);
                obj.amplitude = ones(1, samples);
                obj.phase = zeros(1, samples);
            elseif dim == 2
                samples = round(fieldLength / resolution);
                obj.amplitude = ones(samples, samples);
                obj.phase = zeros(samples, samples);
            end
        end

        function complexField = getComplexField(obj)
            % Returns the complex representation of the field
            complexField = obj.amplitude .* exp(1i * obj.phase);
        end

        function obj = setComplexField(obj, field)
            % Set the amplitude and phase properties from a complex field
            obj.amplitude = abs(field);
            obj.phase = angle(field);
        end

        function obj = resetField(obj)
            % Reset the field amplitude to all 1s and phase to all 0s
            obj.amplitude = ones(size(obj.amplitude));
            obj.phase = zeros(size(obj.phase));
        end

        function obj = applyPhaseShift(obj, phaseShift)
            % Apply a global phase shift
            obj.phase = obj.phase + phaseShift;
        end

        function obj = applyAmplitudeMask(obj, mask)
            % Apply an amplitude mask to the field
            if ~isequal(size(obj.amplitude), size(mask))
                error('Size of mask must match the size of the field amplitude.');
            end
            obj.amplitude = obj.amplitude .* mask;
        end

        function obj = applyPhaseRect(obj, rect_size, phase_shift)
            % Create a rectangle shaped phase shift at the center of the
            % field of a specified size
            if obj.dim == 1
                rect_half_samples = round((rect_size / 2) / obj.resolution);
                middle = round(length(obj.phase) / 2);
                start_idx = max(1, middle - rect_half_samples);
                end_idx = min(length(obj.phase), middle + rect_half_samples);
                obj.phase(start_idx:end_idx) = obj.phase(start_idx:end_idx) + phase_shift;
            elseif obj.dim == 2
                if isscalar(rect_size)
                    rect_half_width = round((rect_size / 2) / obj.resolution);
                    rect_half_height = rect_half_width;
                elseif isvector(rect_size)
                    rect_half_width = round((rect_size(1) / 2) / obj.resolution);
                    rect_half_height = round((rect_size(2) / 2) / obj.resolution);
                end

                middle_x = round(size(obj.phase, 2) / 2);
                middle_y = round(size(obj.phase, 1) / 2);

                start_x = max(1, middle_x - rect_half_width);
                end_x = min(size(obj.phase, 2), middle_x + rect_half_width);
                start_y = max(1, middle_y - rect_half_height);
                end_y = min(size(obj.phase, 1), middle_y + rect_half_height);

                obj.phase(start_y:end_y, start_x:end_x) = obj.phase(start_y:end_y, start_x:end_x) + phase_shift;
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function wdfValues = wdf(obj)
            % Calculate the gradient or difference of the field phase
            if obj.dim == 2
                % 2D case: compute gradient of the phase
                [dx, dy] = gradient(obj.phase, obj.resolution);
                wdfValues = sqrt(dx.^2 + dy.^2); % Magnitude of the gradient
            elseif obj.dim == 1
                % 1D case: compute difference of the phase
                wdfValues = diff(obj.phase) / obj.resolution;
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function fftField = fft(obj)
            % Calculate the Fourier tranform of the complex field
            complexField = obj.getComplexField();
            if obj.dim == 1
                fftField = fftshift(fft(complexField));
            elseif obj.dim == 2
                fftField = fftshift(fft2(complexField));
            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function disp(obj, titleName)
            % Display the amplitude and phase of the field
            if nargin < 2
                titleName = "";
            end

            fig = figure('Color', 'black', 'Position', [657 76 486 846]);

            if obj.dim == 1
                samples = length(obj.amplitude);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                subplot(2, 1, 1);
                plot(xAxis, obj.amplitude, 'LineWidth', 1.5);
                title('Field Amplitude');
                xlabel('Position (m)');
                ylabel('Amplitude');
                grid on;

                subplot(2, 1, 2);
                plot(xAxis, unwrap(obj.phase), 'LineWidth', 1.5);
                title('Field Phase');
                xlabel('Position (m)');
                ylabel('Phase (radians)');
                grid on;

            elseif obj.dim == 2
                samples = size(obj.amplitude, 1);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);
                yAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                subplot(2, 1, 1);
                imagesc(xAxis, yAxis, obj.amplitude);
                colormap(gca, 'gray');
                colorbar('Color', 'white');
                title('Field Amplitude');
                xlabel('x (m)');
                ylabel('y (m)');
                axis equal;
                axis tight;

                subplot(2, 1, 2);
                imagesc(xAxis, yAxis, obj.phase);
                colorbar('Color', 'white');
                colormap(gca, 'default');
                title('Field Phase');
                xlabel('x (m)');
                ylabel('y (m)');
                axis equal;
                axis tight;

            else
                error('Dimensionality must be either 1 or 2.');
            end

            if titleName ~= ""
                sgtitle(titleName, 'FontSize', 18, 'FontWeight', 'bold');
            end

            % Find all axes in the figure and change font color
            ax = findall(fig, 'Type', 'axes');
            for i = 1:length(ax)
                ax(i).XColor = 'white';  % Change X-axis color
                ax(i).YColor = 'white';  % Change Y-axis color
                ax(i).ZColor = 'white';  % Change Z-axis color (if 3D)
                ax(i).Color = 'black';   % Set axes background to black
            end

            % Change title, labels, and other text objects
            textObjs = findall(fig, 'Type', 'text');
            for i = 1:length(textObjs)
                textObjs(i).Color = 'white';  % Change text color
            end
        end

        function dispWDF(obj)
            % Displays the WDF (gradient or difference) of the field phase
            wdfValues = obj.wdf();

            if obj.dim == 1
                samples = length(wdfValues);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                figure('Color', 'white', 'Position', [657 76 486 846]);
                plot(xAxis, wdfValues, 'LineWidth', 1.5);
                title('WDF of Field Phase');
                xlabel('Position (m)');
                ylabel('WDF Value');
                grid on;

            elseif obj.dim == 2
                samples = size(wdfValues, 1);
                xAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);
                yAxis = linspace(-obj.fieldLength / 2, obj.fieldLength / 2, samples);

                figure('Color', 'white', 'Position', [657 76 486 846]);
                imagesc(xAxis, yAxis, wdfValues);
                colormap(gca, 'default');
                colorbar;
                title('WDF of Field Phase');
                xlabel('x (m)');
                ylabel('y (m)');
                axis equal;
                axis tight;

            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function dispFFT(obj)
            % Display the FFT of the field
            fftField = obj.fft();

            if obj.dim == 1
                samples = length(fftField);
                spatialFreq = linspace(-1 / (2 * obj.resolution), 1 / (2 * obj.resolution), samples);

                figure('Color', 'black', 'Position', [657 76 486 846]);

                subplot(2, 1, 1);
                plot(spatialFreq, abs(fftField), 'LineWidth', 1.5);
                title('FFT Amplitude');
                xlabel('Spatial Frequency (1/m)');
                ylabel('Amplitude');
                grid on;

                subplot(2, 1, 2);
                plot(spatialFreq, angle(fftField), 'LineWidth', 1.5);
                title('FFT Phase');
                xlabel('Spatial Frequency (1/m)');
                ylabel('Phase (radians)');
                grid on;

            elseif obj.dim == 2
                samples = size(fftField, 1);
                spatialFreq = linspace(-1 / (2 * obj.resolution), 1 / (2 * obj.resolution), samples);

                figure('Color', 'black', 'Position', [657 76 486 846]);

                subplot(2, 1, 1);
                imagesc(spatialFreq, spatialFreq, abs(fftField));
                colormap(gca, 'gray');
                colorbar;
                title('FFT Amplitude');
                xlabel('Spatial Frequency x (1/m)');
                ylabel('Spatial Frequency y (1/m)');
                axis equal;
                axis tight;

                subplot(2, 1, 2);
                imagesc(spatialFreq, spatialFreq, angle(fftField));
                colormap(gca, 'default');
                colorbar;
                title('FFT Phase');
                xlabel('Spatial Frequency x (1/m)');
                ylabel('Spatial Frequency y (1/m)');
                axis equal;
                axis tight;

            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

    end
end
