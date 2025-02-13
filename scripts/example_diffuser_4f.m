clearvars; clc; close all;
import_dir("../osf/lib");
import_dir("../osf/utils");

%%
sim = Sim(.5e-6, 1e-3, 'paddingRatio', 0, 'lambda', 532e-9, 'dim', 2);

f1 = 100e-3;
f2 = 100e-3;
roughness = 30e-6;
coherence_length = 20e-6;
lens_dist = f1 + f2;
diff_dist_ratio = .5;

sim.addLens(0, f1, 'name', 'Lens 1');
sim.addDiffuser(lens_dist*diff_dist_ratio, roughness, coherence_length, 'name', 'Diffuser');
sim.addLens(lens_dist*(1-diff_dist_ratio), f2, 'name', 'Lens 2');

% sim.propToDist(sim.newField(), 0e-3, 'propMethod', 'rs').disp();
% sim.propToDist(sim.newField(), 25e-3, 'propMethod', 'rs').disp();
% sim.propToDist(sim.newField(), 40e-3, 'propMethod', 'rs').disp();
% sim.propToDist(sim.newField(), 55e-3, 'propMethod', 'rs').disp();
% sim.propToDist(sim.newField(), 75e-3, 'propMethod', 'rs').disp();

sim.prop(sim.newField(), 50e-3, 'verbose', true);

sim.propToDist(sim.newField(), 125e-3).disp();
% sim.propToDist(sim.newField(), 150e-3, 'propMethod', 'rs').disp();
% sim.propToElement(sim.newField(), 'Diffuser', 'propMethod', 'rs').disp();

% analyzePhaseCutoff(sim.elements{2}.apply(sim.newField()).getComplexField(), 532e-9, 1e-7)

%%
import_dir("../osf/lib", "unload");
import_dir("../osf/utils", "unload");
