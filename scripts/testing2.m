clearvars; clc; close all;
addpath('../');
%%

import osf.*;
sim = Sim(1e-6, 1e-3, paddingRatio=2, lambda=532e-9);

sim.addPlane(0e-3, name='Source');
sim.addLens(20e-3, 20e-3, name='Lens 1');

sim.addPlane(20e-3, name='Fourier Plane');
filt = sim.newField().addPhase(pi/2, 'a', [.1e-3 .05e-3], [.08e-3 -.03e-3]).addAmplitude(1, 'c', .05e-3, [-.1e-3 .1e-3]);
sim.addFilter(0, filt, name='Phase Filter', operation='add');

sim.addLens(20e-3, 20e-3, name='Lens 2');
sim.addPlane(20e-3, name='Detector');

dash = Dashboard([2,2]);
dash.show(filt, 'a');
dash.show(filt, 'p');
dash.show(filt, 'a.cross');
dash.show(filt, 'p.wdf', wdfLimit=.5);

