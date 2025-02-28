clearvars; clc; close all;

% Define test grid sizes
grid_sizes = 2.^(5:14);
dx = 10e-6;  % Sampling interval in x (meters per pixel)
dy = 10e-6;  % Sampling interval in y (meters per pixel)
lambda = 500e-9; % Wavelength in meters
num_tests = 10; % Number of times to repeat each test for averaging

% Initialize arrays to store average execution times
time_kz1_avg = zeros(size(grid_sizes));  
time_kz2_avg = zeros(size(grid_sizes)); 

% Loop over different grid sizes
for i = 1:length(grid_sizes)
    Nx = grid_sizes(i);
    Ny = grid_sizes(i);

    % Define frequency grid
    fx = linspace(-1/(2*dx), (1/(2*dx)) - 1/(Nx*dx), Nx);
    fy = linspace(-1/(2*dy), (1/(2*dy)) - 1/(Ny*dy), Ny);
    [FX, FY] = meshgrid(fx, fy);

    % Arrays to store execution times for averaging
    time_kz1_trials = zeros(1, num_tests);
    time_kz2_trials = zeros(1, num_tests);

    % Repeat test `num_tests` times
    for t = 1:num_tests
        % Method 1: Direct sqrt computation
        tic;
        k_z1 = 2 * pi * sqrt((1 / lambda^2) - FX.^2 - FY.^2);
        time_kz1_trials(t) = toc;
        clear k_z1;

        % Method 2: Factorized computation
        tic;
        k_z2 = 2 * pi / lambda * sqrt(1 - lambda^2 .* FX.^2 - lambda^2 .* FY.^2);
        time_kz2_trials(t) = toc;
        clear k_z2;
    end

    % Compute average execution time
    time_kz1_avg(i) = mean(time_kz1_trials);
    time_kz2_avg(i) = mean(time_kz2_trials);

    fprintf('Grid Size: %dx%d | Direct sqrt: %.6f s | Factorized: %.6f s\n', ...
        Nx, Ny, time_kz1_avg(i), time_kz2_avg(i));
    fprintf('Difference: %.6f s\n', time_kz1_avg(i) - time_kz2_avg(i));
end

% Plot averaged execution times
figure;
loglog(grid_sizes, time_kz1_avg, '-o', 'LineWidth', 2, 'DisplayName', 'Direct sqrt Computation');
hold on;
loglog(grid_sizes, time_kz2_avg, '-s', 'LineWidth', 2, 'DisplayName', 'Factorized Computation');

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
title('Execution Time Comparison: Direct sqrt vs. Factorized Computation');
legend show;
