import osf.*;

sim = Sim(5e-6, 5e-3);

sim.addSource();
sim.addLens(100e-3, 100e-3);
sim.addPlane(100e-3, name='Fourier Plane');
sim.addLens(100e-3, 100e-3);
sim.addDetector(100e-3, show=false);

object = sim.newField().setAmplitude(0, 'rect', [2e-3 2e-3]);

fields = sim.prop(object, collect=true);
show(fields);
show(sim);
