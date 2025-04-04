clearvars; clc; close all;
addpath('../');
%%

import osf.Sim;
sim = Sim(1e-6, 1e-3, 'paddingRatio', 1, 'lambda', 532e-9, 'dim', 2);

sim.addPlane(0e-3, 'name', 'Source');
sim.addLens(20e-3, 20e-3, 'name', 'Lens 1');

sim.addPlane(20e-3, 'name', 'Fourier Plane');
filt = sim.newField().addAmplitude(1).img('test_img.png', 'range', [-1 0]);
sim.addFilter(0, filt, 'name', 'Phase Filter');

sim.addLens(20e-3, 20e-3, 'name', 'Lens 2');
sim.addPlane(20e-3, 'name', 'Detector');

rect_size = .2e-3;
rect_offset = .3e-3;
object = sim.newField()...
    .addAmplitude(1)...
    .addPhase(0,      'r', rect_size, [-rect_offset, -rect_offset])...
    .addPhase(pi/2,   'r', rect_size, [rect_offset, -rect_offset])...
    .addPhase(pi,     'r', rect_size, [rect_offset, rect_offset])...
    .addPhase(3*pi/2, 'r', rect_size, [-rect_offset, rect_offset]);

[field, all_fields] = sim.prop(object, 'collect', true);
sim.plotFields(all_fields);

