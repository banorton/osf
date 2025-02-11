clearvars; clc; close all;
import_dir("../osf");

%%
% Implementation of a simple 4f system.
res = 1e-6; field_len = 1e-3;
sim = Sim(res, field_len, 'paddingRatio', 2);

% Add lens for 4f system.
distAfterPreviousElement = 20e-3;
focalLength = 20e-3;
sim.addLens(distAfterPreviousElement, focalLength, 'name', 'Lens 1');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addLens(40e-3, 20e-3, 'name', 'Lens 2');

% Create a blank field with the correct dimensions for the simulation.
rect_width = .2e-3;

% Add a rectangular phase object with a phase shift of pi to the field.
field = sim.newField().applyPhaseRect(rect_width, pi);

% Propagate through all elements and then 10mm past the last element.
% The optional verbose argument makes the function plot the field at all
% critical points in the propagation and display the distances and elements
% at each step.
distAfterLastElement = 20e-3;
sim.prop(field, distAfterLastElement, 'verbose', true);

%%
import_dir("../osf", "unload");
