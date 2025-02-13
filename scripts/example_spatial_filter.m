clearvars; clc; close all;
import_dir("../osf/lib");
import_dir("../osf/utils");

%%
% Create 2D simulation with a resolution of 1um and a view size of 1mm by 1mm.
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2);

% Add elements to the simulation.
sim.addLens(20e-3, 10e-3, 'name', 'objective');
sim.addAperture(10e-3, 'name', 'pinhole', 'circ', 6e-6);
sim.addLens(10e-3, -10e-3, 'name', 'collimating lens');

% Create field with a phase rectangle of width 0.2mm and a phase shift of pi.
field = sim.newField().applyPhaseRect(.2e-3, pi);

% Propagate the field through the system
sim.prop(field, 1e-3, 'verbose', true);

%%
import_dir("../osf/lib", "unload");
import_dir("../osf/utils", "unload");
