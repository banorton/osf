clearvars; clc; close all;
addpath('../');
%%

import osf.*;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 1, 'lambda', 532e-9, 'dim', 2);

ratio = 3/15;
f = 10e-3;
offset = round(2*f*ratio *1e3); % round to nearest mm
offset = offset / 1e3; % convert back to m

sim.addPlane(0, 'name', 'Source');
sim.addLens(f, f, 'name', sprintf('f = %.0fmm', f*1e3));
diffuser = sim.addDiffuser(f-offset, 10e-6, 0, 'name', 'Diffuser');
sim.addLens(f+offset, f, 'name', sprintf('f = %.0fmm', f*1e3));
det = sim.addDetector(f, [40 40], 8e-6, [.5e-3 .5e-3]);

show(sim);
sim.print();

field_i = sim.newField().addAmplitude(1);

corrLen = 5e-6;
diffuser.correlationLength = corrLen;
field_f = sim.prop(field_i);

dash = Dashboard([1 3], figPosition=[0 335 1918 397])...
    .show(diffuser.apply(sim.newField()).phase, 'default', title='Diffuser Phase')...
    .show(field_f, 'default', title='Field @ Detector')...
    .show(field_f, 'wdf', overlayDetector=det, wdfLimit=.15, title='WVD @ Detector');

