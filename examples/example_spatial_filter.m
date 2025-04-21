import osf.*;

% EXAMPLE DESCRIPTION
% This is a simulation for a simple spatial filter. The light initially has
% noise and higher spatial frequencies. It is focused by a 10mm lens,
% filtered at the Fourier plane using a 6um aperture, and collimated with a
% 10mm lens. After the aperture, the light has a smooth gaussian
% distribution.

% Create a simulation with a 1mm by 1mm FOV and 3um resolution.
sim = Sim(3e-6, 1e-3);

% Add a source. If no arugment is entered, the default wavelenth is 632nm.
sim.addSource();

% Add the lenses, pinhole, and collimating lens to the system.
sim.addLens(20e-3, 10e-3, name='f = 10mm', circ=.4e-3);
sim.addAperture(10e-3, name='6um Pinhole', circ=6e-6);
sim.addLens(10e-3, 10e-3, name='f = 10mm', circ=.4e-3);

% Add the detector plane 10mm away from the collimating lens.
sim.addPlane(10e-3, name='Detector');

% Create an object that will add some higher spatial frequencies to the
% system that we can filter out with the spatial filter.
object = sim.newField().setAmplitude(0, 'rect', [.3e-3 .3e-3]);

% Propagate that object field through the system. 'verbose' is an optional
% argument that, when true, will show the field after each element in the
% system.
sim.prop(object, verbose=true);
