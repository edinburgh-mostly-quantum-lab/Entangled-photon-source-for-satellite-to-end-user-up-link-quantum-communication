from neumann_rates import secure_key_rates
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import cycle, product
from collections import defaultdict
plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)

def plot_data(file_name):
    powers = []
    losses = []
    skrs = []
    for power in [1,2.5,5,7.5,10]:#[1,2.5,5,7.5,10]:
        file = [file for file in os.listdir(f"data-with-error/{file_name}") if power == float(file.split("-")[0][:-2])][0]
        with open(f"data-with-error/{file_name}/"+file,"r") as f:
            
            loss, skr, _ = [[float(el) for el in line.split(",")] for line in f.readlines()]
            losses += loss
            skrs += skr
            powers += [power]*len(loss)

    
    powers = [powers[i] for i in range(len(powers)) if skrs[i]>=0]
    losses = [losses[i] for i in range(len(losses)) if skrs[i]>=0]
    skrs = [skrs[i] for i in range(len(skrs)) if skrs[i]>=0]
    plot_area = fig.add_axes([0, 0, 1, 1])
    t=plot_area.tricontourf(losses, powers, np.log10(skrs), levels=10, cmap="Blues")
    fig.colorbar(t)
    plot_area.set_ylabel('Power (mW)',fontsize=10)
    plot_area.set_xlabel('Loss (dB)',fontsize=10)
    fig.savefig(
        f'{__file__.split('.')[0]}.png',
        dpi='figure',
        bbox_inches='tight'
    )

if __name__ == '__main__':
    fig = plt.figure(figsize=(3.5, 1.8))
    run_name = "datarun_17_07_24"
    # the following parameters are obtained by matching to the results in file datarun_17_07_24
    params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    plot_data(run_name)
    