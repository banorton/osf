clearvars; clc; close all;
addpath('../');
%%

import osf.*;

sim = Sim(1e-6, 1e-3, paddingRatio=2, lambda=532e-9);
field = sim.newField().addPhase(pi/2, 'a', [.1e-3 .05e-3], [-.2e-3 0]);

show(field, 'p.wdf', wdfLimit=.2);