clearvars; clc; close all;
import osf.*;

s = Sim(1e-6, 1e-3, paddingRatio=1, lambda=100e-9);

field_i = s.newField().addAmplitude(1).addAmplitude(-1, 'circle', .1e-3, [.1e-3 .1e-3]);
s.addPlane(0, 'name', 'Source');
s.addLens(10e-3, 10e-3, 'name', 'L1');

filtField = s.newField().addPhase(1).addAmplitude(1).addAmplitude(-1, 'a', [20e-6 5e-6]);
s.addFilter(10e-3, filtField, 'operation', 'mult', 'name', 'Filter');

s.addLens(10e-3, 10e-3, 'name', 'L2');
det = s.addDetector(10e-3, [40 40], 10e-6, [.5e-3 .5e-3]);

field_f = s.prop(field_i);

show(s);
s.print();

dash = Dashboard([1 3], figPosition=[47 390 1751 381]);
dash.show(field_i, 'a', title='Object Plane');
dash.show(field_f, 'a', title='Detector Plane Field');
dash.show(field_f, 'wdf', wdfLimit=.15, overlayDetector=det);
show(field_f, plotType='wdf', wdfLimit=.15, overlayDetector=det);
