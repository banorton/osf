function S = AngularSpectrumPropagation_2012_12(U, Xres, Yres, dx, dy, dz,lam, zfx,zfy, Evanescence, directAS, showinfo)

% Tomek Kozacki's AS code: q^2 FFT's (size N x N) instead 1 of (size qN x qN) 
% see Kozacki, T. (2008). Opt.Comm., 281(17), 4219–4223. doi:10.1016/j.optcom.2008.05.023
% and also
% Kozacki T, Falaggis K, & Kujawinska M, Appl. Opt., 51(29), 7080–8 (2012)

n0 = 1;

if (showinfo >0)&&(directAS==0), disp('AngularSpectrumPropagation_2012_12 (q^2 FFTs)'); end
if (showinfo >0)&&(directAS==1), disp('AngularSpectrumPropagation_2013_01 (one FFTs)'); end
if n0 ~= 1, error('ktzerror: method can not be used with not zero n0 '); end
[Ny,Nx] = size(U);

if Nx ~= Ny, 
    error(['ktzerror: the input matrix must be squere: Nx = ' num2str(Nx) ' Ny = ' num2str(Ny)]);
end

if Xres ~= Yres, error('KFerror: both resolutions must be equal');end

if directAS == 0
    mone = ones(1,Nx);
    dfx = 1/Nx/Xres; fx = -Nx/2*dfx : dfx : (Nx/2-1)*dfx; k = 2*pi/lam;
    S = zeros(Nx,Nx);
    nx = ones(Nx,1) * (0:Nx-1); ny = (0:Nx-1).' * ones(1,Nx); 
    counter = 1;

    for ly = 0 : zfy-1 % main algorithm loop
        for lx = 0 : zfx-1 % main algorithm loop

            % STSP 1  - object zero padding
            speedterm = exp(-2*pi*1i*nx*lx/zfx/Nx) .* exp(-2*pi*1i*ny*ly/zfy/Nx);
            ftul = fftshift(fft2(  U.*speedterm  ));

            if Evanescence > 0;
                kernel = exp(-k*abs(dz)*imag( sqrt(1 - lam^2*(mone'*((fx+dfx/zfx*lx).^2)+((fx+dfx/zfy*ly).'.^2)*mone))));
                kernel = kernel.*exp(1i*k*dz*real( sqrt(1 - lam^2*(mone'*((fx+dfx/zfx*lx).^2)+((fx+dfx/zfy*ly).'.^2)*mone))));
            else
                kernel = exp(1i*k*dz*real( sqrt(1 - lam^2*(mone'*((fx+dfx/zfx*lx).^2)+((fx+dfx/zfy*ly).'.^2)*mone))));
            end

            offsetKernel = exp(1i*2*pi*(mone'*((fx+dfx/zfx*lx).*dx)+((fx+dfx/zfy*ly).'.*dy)*mone));
            kernel = kernel .* offsetKernel;

            % STSP 2  - shifted frequency spectra assembly
            S = S + ifft2(ifftshift(  kernel.*ftul  )).*conj(speedterm) /zfx/zfy;
            counter = counter + 1;

        end
    end
end 

if directAS==1
    [Nyold Nxold] = size(U);
    aa = (zfx-1)*Nxold/2; bb = (zfy-1)*Nyold/2; 
    U = padarray(U, [bb aa]);
    [Ny,Nx] = size(U);
    % fs = 1/Xres; 
    k = 2*pi/lam; 
    dfx = 1/Nx/Xres; fx = -Nx/2*dfx : dfx :  (Nx/2-1)*dfx; 
    dfy = 1/Ny/Xres; fy = -Ny/2*dfy : dfy :  (Ny/2-1)*dfy; 
    
    % STSP 1  - object zero padding
    % nx = ones(Nx,1) * (0:Nx-1); ny = (0:Nx-1).' * ones(1,Nx) ; 
    
    ftul = fftshift(fft2(  U  ));
    
    [SQFX SQFY] = meshgrid(fx.^2, fy.^2);
    [FX FY] = meshgrid(fx, fy);
    
    if Evanescence > 0;
        kernel = exp(-k*abs(dz)*imag( sqrt(1 - lam^2*(SQFX+SQFY))));
        kernel = kernel.*exp(1i*k*dz*real( sqrt(1 - lam^2*(SQFX+SQFY))));
    else
        kernel = exp(1i*k*dz*real( sqrt(1 - lam^2*(SQFX+SQFY))));
    end

    offsetKernel = exp(1i*2*pi*( FX.*dx + FY.*dy ));
    kernel = kernel .* offsetKernel;

    % STSP 2  - shifted frequency spectra assembly
    S = ifft2(ifftshift(  kernel.*ftul  ));
    S = S(bb+1:bb+Nyold, aa+1:aa+Nxold);
end 

end
