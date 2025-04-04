import osf.*;

% Create a simulation with 5 µm resolution and 1 mm field of view.
sim = Sim(5e-6, 1e-3);

% Add a source with a wavelength of 532 nm.
sim.addSource(532e-9);

% Define slit separation for the double slit in meters.
slit_separation = .1e-3;

% Create a filter with two narrow horizontal slits.
filter = sim.newField()...
    .setAmplitude(0)...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 -slit_separation/2])...
    .setAmplitude(1, 'rect', [1e-3 .02e-3], [0 slit_separation/2]);

% Add the filter to the system at 10 mm with label "Double Slit".
sim.addFilter(10e-3, filter, name='Double Slit');

% Add an observation plane 50 mm after the slits.
sim.addPlane(50e-3);

% Create a new field and propagate it using 200 steps with live
% visualization of the field as it propagates.
object = sim.newField();
sim.propScan(object, 200, live=true);
