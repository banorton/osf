clearvars; clc; close all;
addpath('../');
%%

import osf.Sim;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

sim.addPlane(0e-3, 'name', 'Source');
sim.addGrating(0, [20, 20], 'type', 'phase');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 1');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 2');
sim.addPlane(20e-3, 'name', 'Detector');

obj = sim.newField().addAmplitude(1);
[~, fields] = sim.prop(obj, 'collect', true);

sim.plotFields(fields);
