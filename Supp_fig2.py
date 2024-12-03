# Author: Thomas Jaeken
# This file generates Figure 2 in the supplementary material of the paper "Entangled photon source for satellite-to-end-user up-link quantum communication"
from neumann_rates import secure_key_rates
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

def plot_model(params, power, loss_range,t_delta=0.5e-9,DC_A=100, DC_B=80, t_dead_A=25e-9, t_dead_B=45e-9):
    """This function plots the secure key rate as a function of the fibre distance and free space loss, 
    based on the model in 10.1103/PhysRevA.104.022406, implemented in the neumann_rates.py file."""
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness,Tcc = params
    setup = secure_key_rates(d=4, t_delta=t_delta, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
    losses = loss_range
    akr = []
    qber = []
    qx = []

    for loss in losses:
        akr += [setup.custom_performance(tcc=Tcc, B=Brightness*power,eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)]
        x = [0.5e-9,Brightness*power]
        meas_cc = setup.__coincidences_measured__(x)
        qber += [setup.__coincidences_erroneous__(x, setup.bit_error) / meas_cc]
        qx += [setup.__coincidences_erroneous__(x, setup.phase_error) / meas_cc]
        
    return np.array(akr), np.array(qber), np.array(qx)

if __name__ == '__main__':
    # The first figure compares the secret key rate for two different heralding efficiencies
    fig = plt.figure(figsize=(3.4, 1.6))
    ax = fig.add_subplot()
    axqber = ax.twinx()
    loss_range = np.linspace(30,55,80)
    params = [9.47893858965945, 9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    akr, qber, qx =plot_model(params, power=5, loss_range=loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(loss_range,akr, color="red", label=r"$\eta_{her}$=20\%")
    axqber.plot(loss_range,qber,color="blue") 
    params = [2.47893858965945, 2.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    akr, qber, qx =plot_model(params, power=5, loss_range=loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(loss_range,akr, "--", color="red", label=r"$\eta_{her}$=100\%")
    axqber.plot(loss_range,qber, "--",color="blue")
    ax.set_yscale('log')
    ax.set_xlabel('Loss (dB)',fontsize=10)
    ax.set_ylabel('Secret key rate (bits/s)',fontsize=10,color='red')
    axqber.set_yscale('linear')
    axqber.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    axqber.set_ylabel('error rate',fontsize=10,color='blue')
    ax.legend(fontsize=6, frameon=False)#(0.8,.5))
    axqber.legend(fontsize=6, frameon=False)#(0.8,.5))
    

    # The second figure compares the secret key rate for two different error rates
    fig = plt.figure(figsize=(3.4, 1.6))
    ax = fig.add_subplot()
    axqber = ax.twinx()

    params = [9.47893858965945, 9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    akr, qber, qx =plot_model(params, power=5, loss_range=loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(loss_range,akr, color="red", label=r"$e_{bit}$=3.8\%")
    axqber.plot(loss_range,qber,color="blue")
    old_akr = akr
    params = [9.47893858965945, 9.89101987420694, 0.0, 0.0, 6135831.8248959305, 1e-09]
    akr, qber, qx =plot_model(params, power=5, loss_range=loss_range,t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(loss_range,akr, "--", color="red", label=r"$e_{bit}$=0\%")
    axqber.plot(loss_range,qber, "--",color="blue")

    ax.set_yscale('log')
    ax.legend(fontsize=6, frameon=False)#(0.8,.5))
    ax.set_xlabel('Loss (dB)',fontsize=10)
    ax.set_ylabel('Secret key rate (bits/s)',fontsize=10,color='red')
    axqber.set_yscale('linear')
    axqber.legend(fontsize=6, frameon=False)#(0.8,.5))
    axqber.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    axqber.set_ylabel('error rate',fontsize=10,color='blue')
    plt.tight_layout()
    plt.show()
