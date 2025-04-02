clearvars; clc; close all;
import osf.*;

sim = Sim(1e-6, 1e-3);

sim.addSource(wavelength=532e-9);
sim.addGrating(0, [10 0]);
sim.addLens(20e-3, 20e-3);
sim.addPlane(20e-3);
sim.addLens(20e-3, 20e-3);
sim.addDetector(20e-3);

show(sim);

rect_width = .2e-3;
% field = sim.newField().addAmplitude(1, 'rect', [.2e-3 .2e-3], [-.3e-3, .3e-3]);
field = sim.newField().addAmplitude(1, 'rect', [.5e-3 .5e-3]);

fields = sim.prop(field, collect=true);
show(fields);
