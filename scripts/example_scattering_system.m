clearvars; clc; close all;
import_dir("../osf/lib");
import_dir("../osf/utils");

%%
sim = Sim(1e-6, 1e-3, 'paddingRatio', 4, 'lambda', 532e-9, 'dim', 2);

roughness = 30e-6; coherence_length = 20e-6;
sim.addDiffuser(0, roughness, coherence_length, 'name', 'diffuser');
sim.addLens(0, 10e-3, 'name', 'Lens 1');

field = sim.newField();
sim.prop(field, 10e-3, 'propMethod', 'rs', 'verbose', true);

%%
import_dir("../osf/lib", "unload");
import_dir("../osf/utils", "unload");
