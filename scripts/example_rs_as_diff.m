clearvars; clc; close all;
import_dir("../osf/lib");
import_dir("../osf/utils");

%%
% Implementation of a simple 4f system.
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2);
sim.addLens(10e-3, 10e-3, 'name', 'Lens 1');
sim.addLens(40e-3, 10e-3, 'name', 'Lens 2');

% Create a blank field with the correct dimensions for the simulation.
field = sim.newField().applyPhaseRect(rect_width, pi);
field.disp("Initial Field");

% Propagate using angular spectrum and Rayleigh-Sommerfeld methods and display final fields.
as_prop = sim.prop(field, 10e-3, 'propMethod', 'as').disp("Using AS Propagation");
rs_prop = sim.prop(field, 10e-3, 'propMethod', 'rs').disp("Using RS Propagation");

%%
import_dir("../osf/lib", "unload");
import_dir("../osf/utils", "unload");
