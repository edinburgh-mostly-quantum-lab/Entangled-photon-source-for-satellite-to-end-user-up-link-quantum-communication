# Author: Thomas Jaeken
# This file generates Figure 5b in "Entangled photon source for satellite-to-end-user up-link quantum communication"
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)
custom_cmap = LinearSegmentedColormap.from_list("red_blue", ["blue", "red"])
fig = plt.figure()
ax = plt.axes(projection='3d')
maxes = []
our_trace = []
files = os.listdir("data/datarun_rerun_tcc_processed")
file_order = np.argsort([float(f[:4]) for f in files])
files = [files[i] for i in file_order]
for i,file in enumerate(files):
    with open("data/datarun_rerun_tcc_processed/"+file,"r") as f:
        data = [[el for el in line.split(",")] for line in f.readlines()]
        headers = data.pop(0)
        data = np.array([[np.nan if el == 'None' else float(el) for el in line ] for line in data])

        rates = data[:,headers.index("Secure Key Rate (bps)")]
        tcc = data[:,headers.index("Coincidence Window (ps)")]
        losses = data[:,headers.index("Loss (dB)")]
        ax.plot3D(losses[::4],tcc[::4],rates[::4],c=custom_cmap(i/len(files)))
        maxes += [[losses[0],tcc[np.argmax(rates)], np.max(rates)]]

        our_trace += [[losses[0],400, rates[np.where(tcc==400)[0][0]]]]
ax.set_ylabel(r"$t_{CC}$ (ps)")
ax.set_xlabel("Loss (dB)")
ax.set_zlabel("Secret key rate (bits/s)")
# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
# ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
maxes = np.array(maxes)
# maxes = maxes[np.argsort(maxes[:,0])]
ax.scatter(*maxes.T,c="green",depthshade=False)
# ax.plot3D(*maxes.T,c="black")
our_trace = np.array(our_trace)
# our_trace = our_trace[np.argsort(our_trace[:,0])]
# ax.plot3D(*our_trace.T,c="grey")
ax.scatter(*our_trace.T,c="black",depthshade=False)
ax.set_box_aspect(aspect=None, zoom=0.8)
ax.view_init(30, -27)
plt.savefig(
    f'{__file__.split('.')[0]}.png',
    dpi='figure',
    bbox_inches='tight'
)
plt.show()