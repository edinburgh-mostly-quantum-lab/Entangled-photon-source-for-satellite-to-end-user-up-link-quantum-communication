
# Author: Thomas Jaeken
# This file generates Figure 2 in the supplementary material of the paper "Entangled photon source for satellite-to-end-user up-link quantum communication"
from neumann_rates import secure_key_rates
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)

def h(x):
    return -x*np.log2(x)-(1-x)*np.log2(1-x)

def plot_model(params, powers, loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9):
    """This function plots the secure key rate as a function of the fibre distance and free space loss, 
    based on the model in 10.1103/PhysRevA.104.022406, implemented in the neumann_rates.py file."""
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc = params
    # Brightness, t_delta,intrinsic_heralding_1550, intrinsic_heralding_780,bit_err, phase_err= params
    setup = secure_key_rates(d=4, t_delta=t_delta, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
    losses = np.linspace(np.min(loss_range),np.max(loss_range),200)

    for power in powers:
        if  not power in [2.5]: continue
        akr = []
        for loss in losses:
            akr += [setup.custom_performance(tcc=Tcc, B=Brightness*power,eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)]

        if DC_B>100:
            ax.plot(losses, akr, '--', c="lightblue")
        else:
            ax.plot(losses, akr, c="lightblue")


def plot_data(file_name, DC_label, correction):
    """This function flicks through the data folder and plots all points by power and loss.
    """
    rates = defaultdict(list)
    undererror = defaultdict(list)
    overerror = defaultdict(list)
    losses = defaultdict(list)
    for file in os.listdir(f"data/{file_name}_processed"):
        if file[0]==".":continue
        with open(f"data/{file_name}_processed/"+file,"r") as f:
            data = [[el for el in line.split(",")] for line in f.readlines()]
            headers = data.pop(0)
            headers[-1] = headers[-1][:-1]  #cut off the \n
            data = np.array([[np.nan if el == 'None' else float(el) for el in line ] for line in data])
            counts = {headers[i]: sum(data[:,i])/sum(data[:,headers.index("Acquisition Time (s)")]) for i in range(headers.index("Acquisition Time (s)")+1,len(headers))}
            power = data[0][headers.index('Power (mW)')]
            losses[power] += [data[0][headers.index('Loss (dB)')]-correction]

            total_CC   = sum(counts[key] for key in ['H1550H780','H1550V780','H1550D780','H1550A780','V1550H780','V1550V780','V1550D780','V1550A780','D1550H780','D1550V780','D1550D780','D1550A780','A1550H780','A1550V780','A1550D780','A1550A780'])

            # The qber and qx will be overestimated at 30db
            qber = (counts['V1550H780']+counts['H1550V780'])/(counts['V1550H780']+counts['H1550V780']+counts['H1550H780']+counts['V1550V780'])
            qx = (counts['A1550D780']+counts['D1550A780'])/(counts['A1550D780']+counts['D1550A780']+counts['D1550D780']+counts['A1550A780'])
            rates[power] += [total_CC/2*(1-h(qber)-h(qx))]

            # The underestimated rate for the errorbar
            qber_under = (counts['V1550H780']+counts['H1550V780']+np.sqrt(counts['V1550H780'])+np.sqrt(counts['H1550V780']))/(counts['V1550H780']+counts['H1550V780']+counts['H1550H780']+counts['V1550V780']+np.sqrt(counts['V1550H780'])+np.sqrt(counts['H1550V780'])-np.sqrt(counts['H1550H780'])-np.sqrt(counts['V1550V780']))
            qx_under = (counts['A1550D780']+counts['D1550A780']+np.sqrt(counts['A1550D780'])+np.sqrt(counts['D1550A780']))/(counts['A1550D780']+counts['D1550A780']+counts['D1550D780']+counts['A1550A780']+np.sqrt(counts['A1550D780'])+np.sqrt(counts['D1550A780'])-np.sqrt(counts['D1550D780'])-np.sqrt(counts['A1550A780']))
            undererror[power] += [(total_CC-np.sqrt(total_CC))/2*(1-h(qber_under)-h(qx_under))]

            # The overestimated rate for the errorbar
            qber_over = (counts['V1550H780']+counts['H1550V780']-np.sqrt(counts['V1550H780'])-np.sqrt(counts['H1550V780']))/(counts['V1550H780']+counts['H1550V780']+counts['H1550H780']+counts['V1550V780']-np.sqrt(counts['V1550H780'])-np.sqrt(counts['H1550V780'])+np.sqrt(counts['H1550H780'])+np.sqrt(counts['V1550V780']))
            qx_over = (counts['A1550D780']+counts['D1550A780']-np.sqrt(counts['A1550D780'])-np.sqrt(counts['D1550A780']))/(counts['A1550D780']+counts['D1550A780']+counts['D1550D780']+counts['A1550A780']-np.sqrt(counts['A1550D780'])-np.sqrt(counts['D1550A780'])+np.sqrt(counts['D1550D780'])+np.sqrt(counts['A1550A780']))
            overerror[power] += [(total_CC+np.sqrt(total_CC))/2*(1-h(qber_over)-h(qx_over))]

    loss_range = [100,0]  ## swapped to force update
    color_gen_data = iter(plt.cm.Blues(np.linspace(1, 0, 5)))
    next(color_gen_data)
    next(color_gen_data)
    for power in np.sort(list(rates.keys())):
        if not power in [5]: continue
        rates[power] = np.array(rates[power])[np.argsort(losses[power])]
        losses[power] = np.array(losses[power])[np.argsort(losses[power])]
        color=next(color_gen_data)
        # ax.errorbar(losses[power], rates[power], yerr=[np.clip(rates[power]-undererror[power],0,None),np.clip(overerror[power]-rates[power],0,None)],linestyle='',marker='o', label=f"{power} mW, {DC_label} DC/s", color=color)
        if DC_label>100:
            ax.scatter(losses[power], rates[power],s=10,label=f"$\sim$ {DC_label} DC/s", color=color,marker="^")
        
        else:
            ax.scatter(losses[power], rates[power],s=10,label=f"$\sim$ {DC_label} DC/s", color=color)
        loss_range[0] = min(losses[power][0],loss_range[0])
        loss_range[1] = max(losses[power][-1],loss_range[1])
    return np.sort(list(rates.keys())), loss_range

if __name__ == '__main__':
    fig = plt.figure(figsize=(3.5, 1.8))
    ax = fig.add_axes([0, 0, 1, 1])
    run_name = 'datarun_17_07_24'
    # the following parameters are the result of matching to datarun_17_07_24
    # params = [1.06439023e+00, 7.24172674e-01, 7.13564918e+00, 7.92265295e+00, 5.00000000e-02, 6.00000000e-02]
    params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 10135831.8248959305, 1e-09]
    powers, loss_range = plot_data(run_name, DC_label=70, correction=0)
    plot_model(params, powers, loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70)

    # increased DC
    run_name = 'datarun_18_07_24_highDC'
    # the following parameters are the result of matching to datarun_18_07_24_highDC
    # params = [0.86439023e+00, 7.24172674e-01, 8.13564918e+00, 7.92265295e+00, 5.00000000e-02, 6.00000000e-02]
    params = [11.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 10135831.8248959305, 1e-09]
    powers, loss_range = plot_data(run_name, DC_label=700, correction=0)
    plot_model(params, powers, loss_range,t_delta=0.4e-9,DC_A=200, DC_B=700, t_dead_A=25e-9, t_dead_B=45e-9)

    ax.set_yscale('log')
    ax.legend(fontsize=10, loc=(0.015,0.03), frameon=False)
    ax.set_xlabel('Loss (dB)',fontsize=10)
    ax.set_ylabel('SKR (bits/s)',fontsize=10)
    plt.show()

