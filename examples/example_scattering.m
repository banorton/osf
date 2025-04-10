import osf.*;

% EXAMPLE DESCRIPTION
% This simulation models light propagation through a simple optical system
% with a diffuser. The setup includes a source, a focusing lens, and a
% random phase diffuser. The diffuser scatters the beam to simulate surface
% roughness. The output is observed at a detector plane.

% Initialize simulation with 5µm resolution and 1.8mm field of view.
sim = Sim(5e-6, 1.8e-3);

% Add a default coherent light source at the start of the system.
sim.addSource();

% Add a lens with 15mm focal length located at position 0.
sim.addLens(0, 15e-3);

% Define diffuser surface properties: roughness and correlation length.
roughness = 20e-6;
correlationLength = 10e-6;

% Add a phase diffuser 10mm after the source.
sim.addDiffuser(10e-3, roughness, correlationLength);

% Add a detector plane 10mm after the diffuser.
sim.addPlane(20e-3, name='Detector Plane');

% Propagate a new field through the optical system.
field = sim.prop(sim.newField());

% Display the intensity at the detector plane.
show(field.intensity, title='Detector');
