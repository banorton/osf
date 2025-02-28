clearvars; clc; close all;

% Define test grid sizes
grid_sizes = 2.^(5:14); % Grid sizes from 2^8 (256) to 2^13 (8192)
num_tests = 10; % Number of times to repeat each test for averaging
dx = 10e-6;  % Sampling interval in x (meters per pixel)
dy = 10e-6;  % Sampling interval in y (meters per pixel)

% Initialize arrays to store average execution times
time_pre_square_avg = zeros(size(grid_sizes));  
time_square_after_avg = zeros(size(grid_sizes)); 

% Loop over different grid sizes
for i = 1:length(grid_sizes)
    Nx = grid_sizes(i);
    Ny = grid_sizes(i);
    
    % Define spatial frequency vectors
    fx = linspace(-1/(2*dx), (1/(2*dx)) - 1/(Nx*dx), Nx);
    fy = linspace(-1/(2*dy), (1/(2*dy)) - 1/(Ny*dy), Ny);

    % Arrays to store execution times for averaging
    time_pre_square_trials = zeros(1, num_tests);
    time_square_after_trials = zeros(1, num_tests);

    % Repeat test `num_tests` times
    for t = 1:num_tests
        % Method 1: Pre-square before meshgrid
        tic;
        [freqX2, freqY2] = meshgrid(fx.^2, fy.^2);
        [freqX, freqY] = meshgrid(fx, fy);
        time_pre_square_trials(t) = toc;
        clear freqX2 freqY2 freqX freqY;

        % Method 2: Squaring After meshgrid
        tic;
        [freqX, freqY] = meshgrid(fx, fy);
        freqX2 = freqX.^2;
        freqY2 = freqY.^2;
        time_square_after_trials(t) = toc;
        clear freqX2 freqY2 freqX freqY;
    end

    % Compute average execution time
    time_pre_square_avg(i) = mean(time_pre_square_trials);
    time_square_after_avg(i) = mean(time_square_after_trials);

    fprintf('Grid Size: %dx%d | Pre-square: %.6f s | Square after: %.6f s\n', ...
        Nx, Ny, time_pre_square_avg(i), time_square_after_avg(i));
    fprintf('Difference: %.6f s\n', time_pre_square_avg(i) - time_square_after_avg(i));
end

% Plot averaged execution times
figure;
loglog(grid_sizes, time_pre_square_avg, '-o', 'LineWidth', 2, 'DisplayName', 'Pre-square before meshgrid');
hold on;
loglog(grid_sizes, time_square_after_avg, '-s', 'LineWidth', 2, 'DisplayName', 'Square after meshgrid');

% Set base-2 tick marks and labels
set(gca, 'XTick', grid_sizes);  % Ensure ticks are at powers of 2
set(gca, 'XTickLabel', strcat('2^{', string(log2(grid_sizes)), '}')); % Convert to 2^power format
set(gca, 'XScale', 'log', 'YScale', 'log'); % Log-log scale

% Adjust grid to align with base-2 spacing
grid on;
ax = gca;
ax.XMinorGrid = 'off'; % Disable minor grid lines
ax.YMinorGrid = 'off'; % Disable minor grid lines

xlabel('Grid Size (Nx = Ny)');
ylabel('Average Execution Time (seconds)');
title('Execution Time Comparison: Pre-square vs. Square After');
legend show;
