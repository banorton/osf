clearvars; clc; close all;
import_dir("../osf/lib");

%%

sim = Sim(.5e-6, 1e-3, 'paddingRatio', 1.5, 'lambda', 532e-9, 'dim', 2);

sim.addPlane(0, 'name', '');
sim.addLens(50e-3, 50e-3, 'name', 'f = 50mm');
diffuser = sim.addDiffuser(43e-3, 1e-6, 10e-6, 'name', 'Diffuser');
sim.addLens(57e-3, 50e-3, 'name', 'f = 50mm');
sim.addPlane(50e-3, 'name', 'Detector');

sim.print();
sim.disp();

field = sim.newField().addAmplitudeRect(.2e-3, 10);
field.disp('Inital Field');

% for i = 5e-6:5e-6:5e-6
%     diffuser.roughness = i;
%     diffuser.correlationLength = i;
%     diffuser.print();
%     titleStr = sprintf('roughness: %.1fmm, coherenceLength: %.1fmm', diffuser.roughness*1e6, diffuser.correlationLength*1e6);
%     sim.prop(field).disp(titleStr);
%     pause(.001);
% end

field_out = sim.prop(field, 'verbose', true);
field_out.cross();

%%

import_dir("../osf/lib", "unload");
