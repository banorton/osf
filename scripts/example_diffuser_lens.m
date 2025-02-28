clearvars; clc; close all;
addpath('../');
%%

% Here we create a simulation with a 1um resolution and a 1mm by 1mm
% window. It has a padding ratio of 2 meaning that, when propagation, it
% will add 0 padding 2 times the size of length/resolution on all sides.
% The wavelength is green and it is in 2 dimenions.
import osf.Sim;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

% The properties of the diffuser are defined and it is added to the system
% 0mm away from the origin.
roughness = 30e-6; correlationLength = 20e-6;
sim.addDiffuser(0, roughness, correlationLength, 'name', 'Diffuser');

% A lens is also added at 0mm from the origin, but the diffuser will be
% applied to the propagating field first and then the lens.
sim.addLens(0, 10e-3, 'name', 'Lens 1');

% Add a detector plane 10mm from the lens so we can see what it looks like
% after propagating for a bit.
sim.addPlane(10e-3, 'name', 'Detector');

% Show a ray optics diagram of the system.
sim.disp();

% Print the simulation parameters in the console.
sim.print();

% The Sim class has a function called newField which will automatically
% generate a field with the same dimensions and properties as the
% simulation so you don't have to do it manually. Then you can directly
% send the field into the simulation.
field = sim.newField();

% Send the fireld through the system. The optional 'verbose' setting
% causes the simulation to display the amplitude and phase of the field
% after every element in the system.
field = sim.prop(field, 'verbose', true);

