import refnx
from refnx import dataset
import numpy as np
import glob
import matplotlib.pyplot as plt
import models.two_layer as tl
import toolbox as tb
import h5py
import dynesty
import sys

var_list = sys.argv[2].split('_')
ml_number = sys.argv[1]

output_dir = './results/'

file_name = 'ml_{}^'.format(ml_number)
for i in var_list:
    file_name += '{}_'.format(i)
file_name = file_name[:-1]

data_dir = './data/ml_{}'.format(ml_number)
data_files = sorted(glob.glob('{}/*_3c.dat'.format(data_dir)))
samples = len(data_files)
datasets = []
for i in data_files:
    d = dataset.ReflectDataset(i)
    datasets.append(dataset.ReflectDataset([d.x, d.y, d.y_err]))

head = {"C": 10, "H": 18, "O": 8, "N": 1, "P": 1}
tail = {"C": 15 * 2, "H": 15 * 4 + 2}

b_head = tb.get_scattering_length(head, 12.5)
b_tail = tb.get_scattering_length(tail, 12.5)
lipids = []
for i in range(samples):
    lipids.append(
        tl.TwoLayer([b_head, b_tail], name='sample{}'.format(i+1)))

air = refnx.reflect.SLD(0, "air")
water = refnx.reflect.SLD(9.45, "h2o")
structures = []
for i in range(samples):
    structures.append(air(0, 0) | lipids[i] | water(0, 3.3))

lipids[0].thick_t.setp(18., vary=True, bounds=(14, 20))
lipids[0].thick_h.setp(10., vary=False)
lipids[0].mol_vol_h.setp(319, vary=False)
lipids[0].mol_vol_t.setp(829, vary=False)
structures[0][-1].rough.setp(3., vary=True, bounds=(2.9, 5))
lipids[0].rough_h_t.constraint = structures[0][-1].rough
lipids[0].rough_t_a.constraint = structures[0][-1].rough
lipids[0].phi_t.setp(0., vary=True, bounds=(0., 0.25))
lipids[0].phi_h.setp(0.5, vary=True, bounds=(0.5, 0.95))

lipids[0].solv_sld.constrain = structures[0][-1].sld.real

for j in range(1, samples):
    lipids[j].thick_h.constraint = lipids[0].thick_h
    if "tt".format(j) in var_list:
        lipids[j].thick_t.setp(18., vary=True, bounds=(14, 20))
    else:
        lipids[j].thick_t.constraint = lipids[0].thick_t
    if "phih".format(j) in var_list:
        lipids[j].phi_h.setp(0.5, vary=True, bounds=(0.2, 0.95))
    else:
        lipids[j].phi_h.constraint = lipids[0].phi_h
    if "rough".format(j) in var_list:
        structures[j][-1].rough.setp(3., vary=True, bounds=(2.9, 5))
    else:
        structures[j][-1].rough.constraint = structures[0][-1].rough
    if "phit".format(j) in var_list:
        lipids[j].phi_t.setp(0., vary=True, bounds=(0., 0.25))
    else:
        lipids[j].phi_t.constaint = lipids[0].phi_t
    lipids[j].rough_h_t.constraint = structures[j][-1].rough
    lipids[j].rough_t_a.constraint = structures[j][-1].rough
    lipids[j].solv_sld.constraint = structures[0][-1].sld.real
    lipids[j].mol_vol_h.constraint = lipids[0].mol_vol_h
    lipids[j].mol_vol_t.constraint = lipids[0].mol_vol_t

models = []

for i in range(samples):
    models.append(refnx.reflect.ReflectModel(structures[i]))

models[0].bkg.setp(datasets[0].y.min(), vary=False)
models[0].scale.setp(1, vary=False)
for i in range(1, samples):
    models[i].scale.constraint = models[0].scale
    models[i].bkg.setp(datasets[i].y.min(), vary=False)


objectives = []

for i in range(samples):
    objectives.append(refnx.analysis.Objective(
        models[i], datasets[i], 
        transform=refnx.analysis.Transform("YX4")))

global_objective = refnx.analysis.GlobalObjective(objectives)

fitter = refnx.analysis.CurveFitter(global_objective)
fitter.fit('differential_evolution', polish=False, seed=1)
result = fitter.fit('L-BFGS-B')

det = np.linalg.det(np.linalg.inv(result.hess_inv.todense()))
m = len(global_objective.varying_parameters())
logz = np.log(np.math.factorial(m) * (4 * np.pi) ** m / (np.sqrt(det))) - global_objective.chisqr() / 2
out_file = open(output_dir + file_name + ".txt", "w")
out_file.write('{}'.format(logz))
out_file.close()
print(global_objective)
