import osf.*;

sim = Sim(10e-6, 5e-3);
sim.addSource(532e-9);
sim.addLens(100e-3, 100e-3);
sim.addPlane(100e-3, name='Fourier Plane');
sim.addLens(100e-3, 100e-3);
sim.addDetector(100e-3, show=false);

show(sim);

sim.removeElement("Lens 1");

show(sim);

sim.elements{3}.focalLength = 50e-3;

show(sim);

object = sim.newField().setAmplitude(0, 'rect', [2e-3 2e-3]);
fields = sim.prop(object, collect=true);
