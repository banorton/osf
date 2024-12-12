clearvars; clc; close all;
import_dir("../osf");

%%
% Implementation of a simple 4f system.
res = 1e-6; field_len = 1e-3;
sim = Sim(res, field_len, 'paddingRatio', 2);

% Add lens for 4f system.
sim.addLens(10e-3, 10e-3, 'name', 'Lens 1');
sim.addLens(40e-3, 10e-3, 'name', 'Lens 2');

% Create a blank field with the correct dimensions for the simulation.
rect_width = .2e-3;
field = sim.newField().applyPhaseRect(rect_width, pi);

% Propagate through all elements and then 10mm past the last element.
% The optional verbose argument makes the function plot the field at all
% critical points in the propagation and display the distances and elements
% at each step.
distAfterLastElement = 10e-3;
sim.prop(field, distAfterLastElement, 'verbose', true);

%%
import_dir("../osf", "unload");
