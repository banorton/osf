clearvars; clc; close all;
import osf.*;

%% -----------SYSTEM DESIGN-----------

fieldLength = 1e-3;
spatialResolution = 3e-6;

ra = fieldLength;
rb = 3.1e-6;
A = 1;
thetaSpan = linspace(0, 2*pi, 10);
BSpan = linspace(0.1, 3.5, 10);
[thetaGrid, BGrid] = meshgrid(thetaSpan, BSpan);
phaseMax = 3*pi/2;

sim = Sim(spatialResolution, fieldLength);
sim.addSource(640e-9);
sim.addLens(10e-3, 10e-3, circ=.5e-3);
sim.addPlane(10e-3, name='FP');
sim.addFilter(0, sim.newField(), operation='add', name='P');
sim.addFilter(0, sim.newField(), operation='mult', name='A');
sim.addLens(10e-3, 10e-3, circ=.5e-3, name='L2');
sim.addPlane(10e-3, name='Detector');
pf = sim.elements{4};
af = sim.elements{5};
show(sim);

baseField = sim.newField();
inputField = sim.newField();
inputField.phase = rand(baseField.size) * phaseMax;

visibility = zeros(size(BGrid));
peakIntensity = zeros(size(BGrid));
psi = zeros(size(BGrid));
Cmag = zeros(size(BGrid));

%% --------------PROPAGATION--------------

for idx = 1:numel(BGrid)
    clc;
    fprintf("%d/%d\n", idx, numel(BGrid));
    B = BGrid(idx);
    BTheta = thetaGrid(idx);
    dtheta = pi/2;
    ATheta = BTheta + dtheta;

    C = B/A * exp(1j * dtheta) - 1;
    psi(idx) = angle(C);
    Cmag(idx) = abs(C);

    pf.field = sim.newField()...
        .setAmplitude(0)...
        .setPhase(BTheta, 'c', rb)...
        .setPhase(ATheta, 'a', [rb ra]);

    af.field = sim.newField()...
        .addAmplitude(B, 'c', rb)... 
        .addAmplitude(A, 'a', [rb ra])...
        .addPhase(1);

    output = sim.prop(inputField);
    intensity = abs(output.amplitude).^2;
    roi = osf.utils.genMask(output.size, 'circle', round(sim.samples/2*.6));
    [Imin, Imax] = bounds(intensity(roi));

    % Dashboard([2 3], figPosition=[72 127 1751 760])...
    %     .showOn(1, inputField)...
    %     .showOn(2, output)...
    %     .showOn(4, inputField, 'p')...
    %     .showOn(5, output, 'p')...
    %     .showOn(3, pf.field, 'p');
    % pause(.01);
    % close all;

    visibility(idx) = (Imax - Imin) / (Imax + Imin);
    peakIntensity(idx) = Imax;
end

%% -----------------PLOTS-----------------

psiVals = mod(psi + pi, 2*pi); % wrapped psi from BTheta
CmagVals = abs(BSpan);               % |C| from BGrid

figure('color', 'white');
imagesc(psiVals, CmagVals, visibility);
axis xy;
xlabel('\psi [rad]');
ylabel('|C|');
title('Visibility (\phi_{max} = 3\pi/2)');
colorbar;
colormap gray;

figure('color', 'white');
imagesc(psiVals, CmagVals, peakIntensity);
axis xy;
xlabel('\psi [rad]');
ylabel('|C|');
title('Peak Irradiance (\phi_{max} = 3\pi/2)');
colorbar;
colormap gray;
