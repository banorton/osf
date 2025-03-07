function W = wdf(signal)
    N = length(signal);

    % Initialize Wigner distribution matrix
    W = zeros(N, N);

    % Compute WDF
    for t = 1:N  % Time index
        for tau = -N/2:N/2-1  % Lag index (symmetric around zero)
            n1 = mod(t + tau, N) + 1;  % Circular indexing
            n2 = mod(t - tau, N) + 1;

            % Compute WDF using the signal and its shifted version
            W(t, tau + N/2 + 1) = signal(n1) * conj(signal(n2));
        end
    end

    % Apply Fourier transform along the lag dimension
    W = fftshift(fft(W, [], 2), 2);  % Compute FFT and shift frequency axis

    % Take real part (WDF should be real-valued)
    W = real(W);
end
