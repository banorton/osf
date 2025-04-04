import osf.*;

% EXAMPLE DESCRIPTION
% This is an example showing how to build a very simple 4f system and
% propagate a simple rectangular object through the system.

% Start by creating a simulation and defining the spatial resolution and
% size of the window respectively.
sim = Sim(10e-6, 5e-3);

% Add the source with the wavelength in meters. If nothing is entered, it
% will default to 632nm. The source plane will always be at the origin.
% When an object is propagated through the system, the source illumination
% will be applied to the input object/field, and propagate starting at the
% source plane.
sim.addSource(532e-9);

% When adding an element to the simulation, the distance from the last
% element always comes first, 100e-3. The Lens element requires a second
% input which is the focal length of the lens.
sim.addLens(100e-3, 100e-3);

% A plane at the Fourier plane which doesn't affect the propagating
% wave in any way. Planes can be used for observation. I give it a name
% which will show up if the system is plotted or if the wave at different
% propagation planes is shown.
sim.addPlane(100e-3, name='Fourier Plane');

% This is the collimating lens.
sim.addLens(100e-3, 100e-3);

% Again, adding elements in this fashion always requires that the first
% input is the distance from the last element, 100mm in this case. A camera
% resolution and pixel pitch are optional 2nd and 3rd inputs. The optional
% argument show is used to tell the system not to display what will be seen
% on the detector when the propagating plane reaches the element.
sim.addDetector(100e-3, show=false);

% Now a simple object is created by first creating a new field. Using
% .newField() on a simulation object creates a plane wave (1 amplitude
% everywhere and 0 phase everywhere) with the same size and resolution as
% the simulation. .setAmplitude() is used on the returned Field object to
% set a 2mm by 2mm rectangle equal to 0 at the center.
object = sim.newField().setAmplitude(0, 'rect', [2e-3 2e-3]);

% To propagate a field through all the elements in a system, the .prop()
% function is used. The first input is the field to be propagated. The
% optional argument collect is used to store and return the propagating
% field after each element.
fields = sim.prop(object, collect=true);

% The show() function can be used on objects like Sim, Field, a collection
% of Fields, or just any matrix. Here it is used to show the field at each
% plane in the system and the layout of the elements in the system.
show(fields);
show(sim);
