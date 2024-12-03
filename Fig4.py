from neumann_rates import secure_key_rates
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import cycle, product
from collections import defaultdict
from functools import partial
plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)

def plot_model(params, powers, loss_range,t_delta=0.5e-9,DC_A=100, DC_B=80, t_dead_A=25e-9, t_dead_B=45e-9):
    """This function plots the secure key rate as a function of the fibre distance and free space loss, 
    based on the model in 10.1103/PhysRevA.104.022406, implemented in the neumann_rates.py file."""
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc = params

    setup = secure_key_rates(d=4, t_delta=t_delta, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
    losses = loss_range
    for power in powers:
        akr = []
        meas_cc = []
        heralding = []
        qber = []
        qx = [] 
        for loss in losses:
            akr += [setup.custom_performance(tcc=Tcc, B=Brightness*power,eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)]
        
        return akr

def plot_data(file_name, losses, power, time):
    """This function flicks through the data folder and plots all points by power and loss."""
    rates = np.zeros_like(time,dtype=float)
    errors = np.zeros_like(time,dtype=float)
    file = f"data-with-error/{file_name}/5.mW-averagedSKR-and-error.csv"
    with open(file,"r") as f:
        loss_range, skr, error = [[float(el) for el in line.split(",")] for line in f.readlines()]
        for i, loss in enumerate(loss_range):
            if np.min(losses)>loss or loss>np.max(losses): continue
            index = np.argmin((losses[:int(len(losses)/2)]-loss)**2)
            rates[index] = skr[i]
            errors[index] = error[i]
        for i, loss in enumerate(loss_range):
            if np.min(losses)>loss or loss>np.max(losses): continue
            index = np.argmin((losses[int(len(losses)/2):]-loss)**2)+int(len(losses)/2)
            rates[index] = skr[i]
            errors[index] = error[i]
    return rates, errors


if __name__ == '__main__':
    fig = plt.figure(figsize=(3.5, 2.6))
    ax = fig.add_subplot()

    # We extracted the loss profile from figure 2 in https://doi.org/10.1038/nature23675 and fitted a function to it.
    time = np.linspace(0,360,360)
    loss_behaviour = lambda x,freq,a,b,c,d: a*np.log10(b**2+c**2-2*b*c*np.cos(freq*(x-d))) - 8.9
    loss_overpass = partial(loss_behaviour,freq=-8.13786355e-05,  a=1.26519905e+01,  b=6.57313691e+03,  c=6.52804180e+03,d=1.82135169e+02)
    axloss = ax.twinx()
    axloss.plot(time,loss_overpass(time),'--', color='black')
    axloss.set_ylabel('Loss (dB)', color='black')

    # 0km
    run_name = 'datarun_02_07_24'
    # the following parameters are obtained by matching to the results in file datarun_02_07_24
    # params = [7.761328918344407, 7.542259886343475, 0.0335677812551474, 0.02531375770896907, 10826895.017621633, 4e-10]
    params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 12135831.8248959305, 4e-010]
    akr=np.array(plot_model(params, powers=[5], loss_range=loss_overpass(time),t_delta=0.3e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9))
    ax.plot(time[np.where(akr>=0)],akr[np.where(akr>=0)],color="darkblue")
    datapoints,yerr = plot_data(run_name,loss_overpass(time), 5, time)
    ax.scatter(time[np.where(datapoints>0)], datapoints[np.where(datapoints>0)], s=4,marker='o',  color="darkblue")

    # 10km
    run_name = 'datarun_17_07_24'
    # the following parameters are obtained by matching to the results in file datarun_17_07_24
    #intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc
    # params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 12135831.8248959305, 1e-09]
    akr=np.array(plot_model(params, powers=[5], loss_range=loss_overpass(time),t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9))
    ax.plot(time[np.where(akr>=0)],akr[np.where(akr>=0)],color="lightblue")
    datapoints,yerr = plot_data(run_name,loss_overpass(time), 5, time)
    ax.scatter(time[np.where(datapoints>0)], datapoints[np.where(datapoints>0)], s=4,marker='o',  color="lightblue")

    ax.set_xlabel('Time (s)',fontsize=10)
    ax.set_ylabel('Secret key rate (bits/s)',fontsize=10)#,color='blue')
    plt.tight_layout()
    plt.show()
    plt.savefig(
        f'{__file__.split('.')[0]}.png',
        dpi='figure',
        bbox_inches='tight'
    )
