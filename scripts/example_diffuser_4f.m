clearvars; clc; close all;
import_dir("../osf/lib");

%%
% Here we create a simulation with a 1um resolution and a 1mm by 1mm
% window. It has a padding ratio of 2 meaning that, when propagation, it
% will add 0 padding 2 times the size of length/resolution on all sides.
% The wavelength is green and it is in 2 dimenions.
sim = Sim(.5e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

% Defining properties of the diffuser.
roughness = 30e-6; coherence_length = 20e-6;

% Add elements to the system.
sim.addPlane(0, 'name', '')
sim.addLens(0, 100e-3, 'name', 'f = 100mm');
sim.addDiffuser(99e-3, roughness, coherence_length, 'name', 'Diffuser');
sim.addLens(101e-3, 100e-3, 'name', 'f = 100mm');
sim.addPlane(25e-3, 'name', 'Detector')

% Show a ray optics diagram of the system.
sim.disp();

% Print the simulation parameters in the console.
sim.print();
sim.elements{3}.print();

% Propagate the field through the system.
sim.prop(sim.newField(), 'verbose', true);

%%
import_dir("../osf/lib", "unload");
