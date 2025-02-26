clearvars; clc; close all;
import_dir("../osf/lib");

%%

sim = Sim(.5e-6, 1e-3, 'paddingRatio', 1.5, 'lambda', 532e-9, 'dim', 2);

ratio = 3/15;
f = 10e-3;
offset = round(2*f*ratio *1e3); % round to nearest mm
offset = offset / 1e3; % convert back to m

sim.addPlane(0, 'name', 'Source');
sim.addLens(f, f, 'name', sprintf('f = %.0fmm', f*1e3));
diffuser = sim.addDiffuser(f-offset, 10e-6, 0, 'name', 'Diffuser');
sim.addLens(f+offset, f, 'name', sprintf('f = %.0fmm', f*1e3));
sim.addPlane(f, 'name', 'Detector');

sim.print();
sim.disp();

field = sim.newField();
field.disp('title', 'Inital Field');

corrLen = 0e-6;
diffuser.correlationLength = corrLen;
diffuser.print();
titleStr = sprintf('roughness: %.1fmm, correlation length: inf', diffuser.roughness*1e6);
sim.prop(field).disp('title', titleStr, 'cross', true);

for corrLen = 30e-6:-5e-6:1e-6
    diffuser.correlationLength = corrLen;
    diffuser.print();

    titleStr = sprintf('roughness: %.1fmm, correlation length: %.1fmm', diffuser.roughness*1e6, diffuser.correlationLength*1e6);

    sim.prop(field).disp('title', titleStr, 'cross', true);
end

%%

import_dir("../osf/lib", "unload");
