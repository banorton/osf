clearvars; clc; close all;
addpath('../');
%%

import osf.Sim;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

filt = sim.newField().addPhase(pi/2, 'c', .05e-3).addAmplitude(0);

sim.addPlane(0e-3, 'name', 'Source');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 1');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addFilter(0, filt, 'name', 'Phase Filter', 'operation', 'add');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 2');
sim.addPlane(20e-3, 'name', 'Detector');

field = sim.newField().addAmplitude(.0001, 'a', [.2e-3 .1e-3], [-.3e-3, .3e-3]);
field = field.addPhase(pi/2, 'a', [.2e-3 .1e-3], [-.3e-3, .3e-3]);

sim.prop(field, 'verbose', true);
