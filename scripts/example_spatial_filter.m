clearvars; clc; close all;
addpath('../');

%%

% Here we create a simulation with a 1um resolution and a 1mm by 1mm
% window. It has a padding ratio of 2 meaning that, when propagation, it
% will add 0 padding 2 times the size of length/resolution on all sides.
% The wavelength is green and it is in 2 dimenions.
import osf.Sim;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

% Add elements to the simulation.
sim.addPlane(0, 'name', 'Source');
sim.addLens(20e-3, 10e-3, 'name', 'f = 10mm');
sim.addPlane(10e-3, 'name', 'Focal Plane');
sim.addAperture(0, 'name', '6um Pinhole', 'circ', 6e-6);
sim.addLens(10e-3, 10e-3, 'name', 'f = 10mm');
sim.addPlane(20e-3, 'name', 'Detector');

% Show a ray optics diagram of the system.
sim.disp();

% Print the simulation parameters in the console.
sim.print();

% Create field with a phase rectangle of width 0.2mm and a phase shift of pi.
field = sim.newField().addPhaseRect(.2e-3, pi);

% Propagate the field through the system.
sim.prop(field, 'verbose', true);

