import osf.*;

% EXAMPLE DESCRIPTION
% Binary amplitude grating with lens and a detector placed at the Fourier
% plane.

% Initialize simulation with 5µm resolution and 1mm x 1mm FOV.
sim = Sim(5e-6, 1e-3);

% Add a default coherent source plane.
sim.addSource();

% Add a binary grating 0mm from the source with 0 lines/mm along x, and 20
% lines/mm along y.
sim.addGrating(0, [0 20], gratingType='binary');

% Add a lens 20mm after the grating with 20mm focal length.
sim.addLens(20e-3, 20e-3);

% Add a detector 20mm after the lens.
sim.addDetector(20e-3);

% Propagate a plane wave through the system.
field = sim.prop(sim.newField());

% Display the x and y cross section of the field amplitude after
% propagation.
show(field, 'cross');
