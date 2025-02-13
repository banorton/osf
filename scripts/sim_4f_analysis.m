clearvars; clc; close all;
import_dir("../osf/lib");
import_dir("../osf/utils");

%%
res = 1e-6; field_len = 1e-3;
sim = Sim(res, field_len, 'paddingRatio', .5);

sim.addLens(20e-3, 20e-3, 'name', 'Lens 1');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addLens(40e-3, 20e-3, 'name', 'Lens 2');

sim.disp();

% field = sim.newField().applyPhaseRect(.2e-3, pi);
% field = sim.newField();

% sim.prop(field, 20e-3, 'verbose', true);
% test = sim.propToElement(sim.newField(), 'Fourier Plane');
% test.disp();
% test.phaseCross();

%%
import_dir("../osf/lib", "unload");
import_dir("../osf/utils", "unload");
