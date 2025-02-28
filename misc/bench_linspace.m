clearvars; clc; close all;

% Define test grid sizes (powers of 2)
grid_sizes = 2.^(5:15); % Grid sizes from 2^5 (32) to 2^14 (16384)
num_tests = 10; % Number of     times to repeat each test for averaging
dx = 10e-6;  % Sampling interval in x (meters per pixel)
dy = 10e-6;  % Sampling interval in y (meters per pixel)

% Initialize arrays to store average execution times
time_index_avg = zeros(size(grid_sizes));  
time_linspace_avg = zeros(size(grid_sizes)); 

% Loop over different grid sizes
for i = 1:length(grid_sizes)
    Nx = grid_sizes(i);
    Ny = grid_sizes(i);

    % Arrays to store execution times for averaging
    time_index_trials = zeros(1, num_tests);
    time_linspace_trials = zeros(1, num_tests);

    % Repeat test `num_tests` times
    for t = 1:num_tests
        % Method 1: Index-based frequency computation
        tic;
        fx1 = (-Nx/2:Nx/2-1) / (Nx * dx);
        fy1 = (-Ny/2:Ny/2-1) / (Ny * dy);
        [FX1, FY1] = meshgrid(fx1, fy1);
        time_index_trials(t) = toc;
        clear FX1 FY1;

        % Method 2: Optimized linspace computation
        tic;
        fx2 = linspace(-1/(2*dx), (1/(2*dx)) - 1/(Nx*dx), Nx);
        fy2 = linspace(-1/(2*dy), (1/(2*dy)) - 1/(Ny*dy), Ny);
        [FX2, FY2] = meshgrid(fx2, fy2);
        time_linspace_trials(t) = toc;
        clear FX2 FY2;
    end

    % Compute average execution time
    time_index_avg(i) = mean(time_index_trials);
    time_linspace_avg(i) = mean(time_linspace_trials);

    fprintf('Grid Size: %dx%d | Index Method: %.6f s | Linspace Method: %.6f s\n', ...
        Nx, Ny, time_index_avg(i), time_linspace_avg(i));
    fprintf('Difference: %.6f s\n', time_index_avg(i) - time_linspace_avg(i));
end

% Plot averaged execution times
figure;
loglog(grid_sizes, time_index_avg, '-o', 'LineWidth', 2, 'DisplayName', 'Index-based Method');
hold on;
loglog(grid_sizes, time_linspace_avg, '-s', 'LineWidth', 2, 'DisplayName', 'Linspace Method');

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
title('Execution Time Comparison: Index-based vs. Linspace Frequency Computation');
legend show;
