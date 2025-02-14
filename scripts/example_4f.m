clearvars; clc; close all;
import_dir("../osf/lib");

%%
% Here we create a simulation with a 1um resolution and a 1mm by 1mm
% window. It has a padding ratio of 2 meaning that, when propagation, it
% will add 0 padding 2 times the size of length/resolution on all sides.
% The wavelength is green and it is in 2 dimenions.
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

% Add elements for 4f system.
distAfterPreviousElement = 20e-3;
focalLength = 20e-3;
sim.addPlane(0e-3, 'name', 'Source');
sim.addLens(distAfterPreviousElement, focalLength, 'name', 'Lens 1');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 2');
sim.addPlane(20e-3, 'name', 'Detector');

% Show ray optics diagram of system.
sim.disp();

% Print the simulation parameters to the console.
sim.print();

% Add a rectangular phase object with a phase shift of pi to the field.
rect_width = .2e-3;
field = sim.newField().applyPhaseRect(rect_width, pi);

% Send the fireld through the system. The optional 'verbose' setting
% causes the simulation to display the amplitude and phase of the field
% after every element in the system.
sim.prop(field, 'verbose', true);

%%
import_dir("../osf/lib", "unload");
