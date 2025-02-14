clearvars; clc; close all;
import_dir("../osf/lib");

%%
res = 1e-6; field_len = 1e-3;
sim = Sim(res, field_len, 'paddingRatio', .5);

sim.addPlane(0, 'name', 'Source');
sim.addLens(40e-3, 20e-3, 'name', 'f = 20mm');
sim.addPlane(20e-3, 'name', 'Fourier Plane');
sim.addLens(20e-3, 20e-3, 'name', 'f = 20mm');
sim.addPlane(40e-3, 'name', 'Detector');

sim.disp();
sim.print();

field = sim.newField();
field = sim.prop(field, 'verbose', true);
% field.phaseCross();

%%
import_dir("../osf/lib", "unload");
