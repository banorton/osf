clearvars; clc; close all;
import osf.*;

s = Sim(1e-6, 1e-3, paddingRatio=1, lambda=100e-9);

% Object
f = s.newField().addAmplitude(1, 'circle', .1e-3, [.1e-3 .1e-3]);

% Source and L1
s.addPlane(0, 'name', 'Source');
s.addLens(10e-3, 10e-3, 'circ', .25e-3, 'name', 'L1');

% Phase Filter
filtField = s.newField().addPhase(1).addAmplitude(1).addAmplitude(-1, 'a', [20e-6 5e-6]);
s.addFilter(10e-3, filtField, 'operation', 'mult', 'name', 'Filter');

% s.addPlane(10e-3, 'name', 'Fourier Plane');

% L2 and Detector
s.addLens(10e-3, 10e-3, 'name', 'L2');
s.addPlane(10e-3, 'name', 'Detector');

filtField.show();
s.prop(f, verbose=true);
% field = s.propToElement(f, 'Fourier Plane');
% show(field, 'a.cross');
