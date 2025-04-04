clearvars; clc; close all;
import osf.*;

sim = Sim(5e-6, 1e-3);

sim.addSource(532e-9);

slit_separation = .1e-3;
filter = sim.newField()...
    .addAmplitude(-1)...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 -slit_separation/2])...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 slit_separation/2]);
sim.addFilter(10e-3, filter, name='Double Slit');
sim.addPlane(50e-3);

object = sim.newField();
sim.propScan(object, 200, live=true);