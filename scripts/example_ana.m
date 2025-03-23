clearvars; clc; close all;
import osf.*;

s = Sim(1e-6, 1e-3, paddingRatio=1, lambda=100e-9);

s.addPlane(0, 'name', 'Source');
% s.addLens(10e-3, 10e-3, 'circ', .25e-3, 'name', 'L1');
s.addLens(10e-3, 10e-3, 'name', 'L1');

filtField = s.newField().addPhase(1).addAmplitude(1).addAmplitude(-1, 'a', [20e-6 5e-6]);
s.addFilter(10e-3, filtField, 'operation', 'mult', 'name', 'Filter');

s.addLens(10e-3, 10e-3, 'name', 'L2');
% s.addPlane(10e-3, 'name', 'Detector');
s.addDetector(10e-3, [30 10], 25e-6, [.5e-3 .5e-3]);

f = s.newField().addAmplitude(1).addAmplitude(-1, 'circle', .1e-3, [.1e-3 .1e-3]);
s.prop(f, verbose=true);
