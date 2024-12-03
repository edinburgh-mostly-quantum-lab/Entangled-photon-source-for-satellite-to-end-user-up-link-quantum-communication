# Author: Thomas Jaeken
# This file generates Figure 5a in "Entangled photon source for satellite-to-end-user up-link quantum communication"
from neumann_rates import secure_key_rates
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from functools import partial
plt.rc('text', usetex=True)
plt.rc('font', family='times', size=10)

def plot_model(params, power, loss_range,t_delta=0.5e-9,DC_A=100, DC_B=80, t_dead_A=25e-9, t_dead_B=45e-9):
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness,Tcc = params
    labels = ['intrinsic_heralding_1550: ', 'intrinsic_heralding_780: ','bit_err: ', 'phase_err: ', 'Brightness: ','t_delta: ']
    print([labels[i]+str(params[i]) for i in range(len(params))]+['DC_A=200, DC_B=70'])
    setup = secure_key_rates(d=4, t_delta=t_delta, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
    losses = loss_range
    akr = []
    opt_akr = []
    opt_power = []
    opt_tcc = []
    for loss in losses:
        akr += [setup.custom_performance(tcc=Tcc, B=Brightness*power,eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)]
        opt = setup.optimize_performance(eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)
        opt_akr += [opt[1]]
        # print(opt_akr)
        opt_power += [opt[0][1]/Brightness]
        opt_tcc += [opt[0][0]]
    return np.array(akr), np.array(opt_akr), np.array(opt_power), np.array(opt_tcc)

if __name__ == '__main__':
    fig = plt.figure(figsize=(3.4, 2.6))
    ax = fig.add_subplot()

    # plot the overpass fitted from '/Users/thomasjaeken/Heriot-Watt University Team Dropbox/RES_EPS_EMQL/projects/Optical ground station/__experimental-work__/__BBM92__/loss profile for micius uplink/loss profile for micius uplink.csv'
    time = np.concatenate((np.linspace(0,30,600),np.linspace(30,140,2000),np.linspace(140,220,2000),np.linspace(220,340,2000),np.linspace(340,360,600)))
    loss_behaviour = lambda x,freq,a,b,c,d: a*np.log10(b**2+c**2-2*b*c*np.cos(freq*(x-d)))-3
    loss_overpass = partial(loss_behaviour,freq=-8.13786355e-05,  a=1.26519905e+01,  b=6.57313691e+03,  c=6.52804180e+03,d=1.82135169e+02)
    axloss = ax.twinx()
    axloss.plot(time,loss_overpass(time),'-', color='black')
    custom_cmap = LinearSegmentedColormap.from_list("purple_green", ["blue", "red"])
    axloss.set_ylabel('Loss (dB)', color='black')
    # the following parameters are obtained by matching to the results in file datarun_17_07_24
    params = [9.47893858965945, 9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    
    akr, opt_akr, opt_power, opt_tcc=plot_model(params, power=5, loss_range=loss_overpass(time),t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(time[np.where(np.array(akr)>=0)],akr[np.where(np.array(akr)>=0)],'--',color=custom_cmap(1-5/max(opt_power[np.where(np.array(opt_akr)>=0)])), label="5mW")
    ax.scatter(time[np.where(np.array(opt_akr)>=0)],opt_akr[np.where(np.array(opt_akr)>=0)],color=custom_cmap(1-opt_power[np.where(np.array(opt_akr)>=0)]/max(opt_power[np.where(np.array(opt_akr)>=0)])), label="P(t)",edgecolor='none',s=3)

    
    ax.legend(fontsize=6, loc=(0,0), frameon=False)#(0.8,.5))
    ax.set_xlabel('Time (s)',fontsize=10)
    ax.set_ylabel('Secret key rate (bits/s)',fontsize=10)#,color='blue')
    plt.tight_layout()
    plt.show()

