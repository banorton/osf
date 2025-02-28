% Get original field dimensions
[Ny_old, Nx_old] = size(U);

% Compute padding amounts
pad_x = (zfx - 1) * Nx_old / 2;
pad_y = (zfy - 1) * Ny_old / 2;

% Apply zero-padding to the input field
U = padarray(U, [pad_y, pad_x]);
[Ny, Nx] = size(U);

% Define wave parameters
k = 2 * pi / lambda;  % Wavenumber

% Compute spatial frequency coordinates
dfx = 1 / (Nx * Xres);
fx = (-Nx/2 : Nx/2 - 1) * dfx;

dfy = 1 / (Ny * Xres);
fy = (-Ny/2 : Ny/2 - 1) * dfy;

% Compute Fourier transform of the padded field
U_freq = fftshift(fft2(U));

% Create spatial frequency square grids
[FX_sqr, FY_sqr] = meshgrid(fx.^2, fy.^2);

% Compute angular spectrum transfer function
if Evanescence > 0
    propagation_kernel = exp(-k * abs(dz) * imag(sqrt(1 - lambda^2 * (FX_sqr + FY_sqr))));
    propagation_kernel = propagation_kernel .* exp(1i * k * dz * real(sqrt(1 - lambda^2 * (FX_sqr + FY_sqr))));
else
    propagation_kernel = exp(1i * k * dz * real(sqrt(1 - lambda^2 * (FX_sqr + FY_sqr))));
end

clear FX_sqr, FY_sqr;

% Create spatial frequency grids
[FX, FY] = meshgrid(fx, fy);

% Apply lateral shift phase correction
shift_phase = exp(1i * 2 * pi * (FX * dx + FY * dy));
propagation_kernel = propagation_kernel .* shift_phase;

clear FX, FY;

% Apply transfer function in frequency domain and transform back
U_propagated = ifft2(ifftshift(propagation_kernel .* U_freq));

% Crop result back to original size
S = U_propagated(pad_y + 1 : pad_y + Ny_old, pad_x + 1 : pad_x + Nx_old);
