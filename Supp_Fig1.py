# Author: Thomas Jaeken
# This file generates Figure 1 in the supplementary material of the paper "Entangled photon source for satellite-to-end-user up-link quantum communication"
from neumann_rates import secure_key_rates
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)

def plot_model(params, fibre_range, loss_range,t_delta=0.5e-9,DC_A=100, DC_B=80, t_dead_A=25e-9, t_dead_B=45e-9, dispersion=0):
    """This function plots the secure key rate as a function of the fibre distance and free space loss, 
    based on the model in 10.1103/PhysRevA.104.022406, implemented in the neumann_rates.py file."""
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc = params

    akr = np.empty((len(fibre_range),len(loss_range)))
    for i,fibre in enumerate(fibre_range):
        for j,loss in enumerate(loss_range):
            setup = secure_key_rates(d=4, t_delta=t_delta+dispersion*fibre, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
            akr[i][j] = setup.custom_performance(tcc=Tcc, B=Brightness,eff_A=intrinsic_heralding_1550+0.2*fibre, eff_B=intrinsic_heralding_780+loss)
    
    X,Y = np.meshgrid(loss_range, fibre_range)
    fig = plt.figure(figsize=(3.5, 2.6))
    ax = fig.add_subplot()
    ax.set_title('Dispersion = '+str(dispersion*1e12)+' ps/km')
    plt.contourf(X, Y, akr, levels=20, cmap='Blues',locator=ticker.LogLocator())
    plt.colorbar()
    ax.set_xlabel('Free space loss (dB)',fontsize=10)
    ax.set_ylabel('Fibre distance (km)',fontsize=10)#,color='blue')
    plt.tight_layout()
    

if __name__ == '__main__':
    fibre_range=np.linspace(0,300,100)
    loss_range=np.linspace(30,50,100)

    # the following parameters are obtained by matching to the results in file datarun_02_07_24
    params = [7.761328918344407, 7.542259886343475, 0.0335677812551474, 0.02531375770896907, 10826895.017621633, 4e-10]
    plot_model(params, fibre_range, loss_range,t_delta=0.5e-9,DC_A=200, DC_B=120, t_dead_A=25e-9, t_dead_B=45e-9)
    plot_model(params, fibre_range, loss_range,t_delta=0.5e-9,DC_A=200, DC_B=120, t_dead_A=25e-9, t_dead_B=45e-9,dispersion=40e-12)
    plt.show()

    
