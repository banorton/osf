import osf.*;

% EXAMPLE DESCRIPTION
% This is an example showing how to build a phase contrast simulation.

% Initialize simulation with 2 µm resolution and 1 mm field of view.
sim = Sim(2e-6, 1e-3);

% Add a source plane labeled "Object Plane".
sim.addSource(name='Object Plane');

% Add a lens with 20mm focal length, 20 mm after the source.
sim.addLens(20e-3, 20e-3);

% Create a phase filter that adds a π/2 phase shift in a small circular
% region at the center.
filt = sim.newField().setAmplitude(0).addPhase(pi/2, 'c', .03e-3);
sim.addFilter(20e-3, filt, operation='add');

% Add a second lens for the 4-f system (20mm focal length, 20mm after
% filter).
sim.addLens(20e-3, 20e-3);

% Add a detector 20mm after the second lens.
sim.addDetector(20e-3);

% Create an object with four rectangular phase zones of different phase
% values.
rect_size = .2e-3;
rect_offset = .2e-3;
object = sim.newField()...
    .addPhase(0,      'r', rect_size, [-rect_offset, -rect_offset])...
    .addPhase(pi/2,   'r', rect_size, [rect_offset, -rect_offset])...
    .addPhase(pi,     'r', rect_size, [rect_offset, rect_offset])...
    .addPhase(3*pi/2, 'r', rect_size, [-rect_offset, rect_offset]);

% Display phase of object and filter.
object.show(title='Object');
filt.show('p', title='Phase Filter');
sim.show();

% Propagate the object field through the system.
sim.prop(object).show();
