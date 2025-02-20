clearvars; clc; close all;
import_dir("../osf/lib");

%%

sim = Sim(.5e-6, 1e-3, 'paddingRatio', 1.5, 'lambda', 532e-9, 'dim', 2);

sim.addPlane(0, 'name', 'Source');
sim.addLens(50e-3, 50e-3, 'name', 'f = 50mm');
diffuser = sim.addDiffuser(43e-3, 10e-6, 5e-6, 'name', 'Diffuser');
sim.addLens(57e-3, 50e-3, 'name', 'f = 50mm');
sim.addPlane(50e-3, 'name', 'Detector');

sim.print();
sim.disp();

field = sim.newField().addAmplitudeRect(.2e-3, 10);
field.disp('title', 'Inital Field');

for corrLen = 0e-6
    diffuser.correlationLength = corrLen;
    diffuser.print();
    % titleStr = sprintf('roughness: %.1fmm, coherenceLength: %.1fmm', diffuser.roughness*1e6, diffuser.correlationLength*1e6);
    titleStr = sprintf('roughness: %.1fmm, coherenceLength: inf', diffuser.roughness*1e6);
    sim.prop(field).disp('title', titleStr, 'cross', true);
end

%%

import_dir("../osf/lib", "unload");
