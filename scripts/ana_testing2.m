clearvars; clc; close all;
addpath('../');
%%

import osf.Sim;
sim = Sim(1e-6, 1e-3, paddingRatio=2, lambda=532e-9);

sim.addPlane(0e-3, name='Source');
sim.addLens(20e-3, 20e-3, name='Lens 1');

sim.addPlane(20e-3, name='Fourier Plane');
filt = sim.newField().addPhase(pi/2, 'c', .05e-3).addAmplitude(0);
sim.addFilter(0, filt, name='Phase Filter', operation='add');

sim.addLens(20e-3, 20e-3, name='Lens 2');
sim.addPlane(20e-3, name='Detector');

field = sim.newField().addAmplitude(1, 'circ', .2e-3);
field = sim.prop(field);


field.disp(cross=true);
[xdata, crossdata] = field.getPhaseCross();
[d, f, t] = wvd(crossdata, 1000000, 'smoothedPseudo', MinThreshold=0, NumFrequencyPoints=1000);
figure; imagesc(t,f,d(1:100,:)); ylabel('Local Spatial Frequency (m^{-1})'); xlabel('x (m)');
