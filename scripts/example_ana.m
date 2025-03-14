clearvars; clc; close all;
import osf.*;

s = Sim(1e-6, 1e-3, paddingRatio=1, lambda=532e-9);

% Object
f = s.newField().addPhase(pi, 'c', .1e-3).addAmplitude(1);

% Source and L1
s.addPlane(0, 'name', 'Source');
s.addLens(10e-3, 10e-3, 'circ', .01e-3, 'name', 'L1');

% Phase Filter
filtField = s.newField().addPhase(pi/2, 'circle', .05e-3);
s.addFilter(10e-3, filtField, 'operation', 'add', 'name', 'Filter');

% L2 and Detector
s.addLens(10e-3, 10e-3);
s.addPlane(10e-3, 'name', 'Detector');

% Propagate field through system.
f_final = s.prop(f, verbose=true);

