variables = ["tt", "rough", "phit", "phih"]

import itertools

VARIABLES = []
for i in range(0, 4):
    for j in itertools.combinations(variables, i + 1):
        VARIABLES.append('_'.join(j))

import glob

ml_folders = glob.glob('data/*')
ml_numbers_analysis = []
for folder in ml_folders:
    if 'ml_' in folder:
        if int(folder.split('_')[1]) not in [34]:
            ml_numbers_analysis.append(int(folder.split('_')[1]))


rule targets:
    input:
        ['results/ml_{}^{}.txt'.format(1, i) for i in VARIABLES],
        ['results/ml_{}^{}.txt'.format(4, i) for i in VARIABLES],
        ['results/ml_{}^{}.txt'.format(12, i) for i in VARIABLES],
        ['results/ml_{}^{}.txt'.format(22, i) for i in VARIABLES],
        ['results/ml_{}^{}.txt'.format(33, i) for i in VARIABLES],
        ['results/ml_{}_chain.h5'.format(i) for i in ml_numbers_analysis],
        ['results/ml_{}_output.txt'.format(i) for i in ml_numbers_analysis]

rule para:
    input:
        'batch.py'
    output:
        'results/ml_{sample}^{seed}.txt'
    run:
        shell("python {input} {wildcards.sample} {wildcards.seed}")

rule analysis:
    input:
        'analysis.py'
    output:
        'results/ml_{number}_chain.h5',
        'results/ml_{number}_output.txt'
    run:
        shell("python {input} {wildcards.number}")
