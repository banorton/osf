clearvars; clc; close all;
import_dir("../osf");

%%
sim = Sim(1e-6, 1e-3, 'paddingRatio', 2, 'lambda', 532e-9, 'dim', 2);

roughness = 30e-6; coherence_length = 20e-6;
sim.addDiffuser(0, roughness, coherence_length, 'name', 'diffuser');

field = sim.newField();

sim.prop(field, 50e-3, 'propMethod', 'rs').disp();

%%
import_dir("../osf", "unload");
