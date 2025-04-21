import osf.*;

% EXAMPLE DESCRIPTION
% This simulation models the classic Young's double-slit experiment using
% Fourier optics. A coherent beam passes through two narrow slits,
% producing an interference pattern at the observation plane. The
% simulation uses a custom filter to define the slits and visualizes the
% field propagation

% Create a simulation with 5µm resolution and 1mm field of view.
sim = Sim(5e-6, 1e-3);

% Add a source with a wavelength of 532nm.
sim.addSource(532e-9);

% Define slit separation for the double slit in meters.
slit_separation = .2e-3;

% Create a filter with two narrow horizontal slits.
filter = sim.newField()...
    .setAmplitude(0)...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 -slit_separation/2])...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 slit_separation/2]);

% Add the filter to the system at 10mm with label "Double Slit".
sim.addFilter(10e-3, filter, name='Double Slit');

% Add an observation plane 50mm after the slits.
sim.addPlane(50e-3, name='Detector Plane');

% Create a new field and propagate it using 200 steps with live
% visualization of the field as it propagates.
object = sim.newField();
sim.propScan(object, 200, live=true);
