import osf.*;

sim = Sim(2e-6, 1e-3);

sim.addSource();
sim.addLens(20e-3, 20e-3);

filt = sim.newField().setAmplitude(0).addPhase(pi/2, 'c', .03e-3);
sim.addFilter(20e-3, filt, operation='add');

sim.addLens(20e-3, 20e-3);
sim.addDetector(20e-3);

rect_size = .2e-3;
rect_offset = .2e-3;
object = sim.newField()...
    .addPhase(0,      'r', rect_size, [-rect_offset, -rect_offset])...
    .addPhase(pi/2,   'r', rect_size, [rect_offset, -rect_offset])...
    .addPhase(pi,     'r', rect_size, [rect_offset, rect_offset])...
    .addPhase(3*pi/2, 'r', rect_size, [-rect_offset, rect_offset]);

object.show();
sim.prop(object);
sim.propScan(object, 200, live=true);