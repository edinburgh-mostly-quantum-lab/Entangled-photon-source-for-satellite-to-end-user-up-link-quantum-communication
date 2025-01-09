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
fig = plt.figure(figsize=(3.5, 2.6))
ax = fig.add_subplot()

def h(x):
    return -x*np.log2(x)-(1-x)*np.log2(1-x)

def find_start_point(file_name):
    """This function flicks through the foder to find the data point with the lowest power and lowest loss to match the model with.
        The other points should follow the model without needing to be matched.
    """
    start_file = ''
    loss = np.inf
    power = np.inf
    for f in os.listdir("data/"+file_name+'_processed'):
        if f[0]==".":continue
        if float(f.split('_')[0]) <= loss and float(f.split('_')[4]) <= power:
            start_file = f
            loss = float(f.split('_')[0])
            power = float(f.split('_')[4])
    return file_name+"_processed/"+start_file, loss, power

def match_model(start_data, loss, power):
    """this functions is meant to find the parameters for the model that fit make it overlap with the first data point.
    It should go through a few iterations of optimization but the first iteration is with hardcoded values.
    ideally this would be a datapoint without loss but it still works if we start at nonzero loss.
    """
    f = open(os.getcwd()+f"/data/{start_data}", 'r')
    data = [[el for el in line.split(",")] for line in f.readlines()]
    headers = data.pop(0)
    headers[-1] = headers[-1][:-1]  #cut off the \n
    data = np.array([[np.nan if el == 'None' else float(el) for el in line ] for line in data])
    Tcc = data[0,headers.index("Coincidence Window (ps)")]*1e-12
    counts = {headers[i]: sum(data[:,i])/sum(data[:,headers.index("Acquisition Time (s)")]) for i in range(headers.index("Acquisition Time (s)")+1,len(headers))}

    total_1550 = sum(counts[key] for key in ['H1550','V1550','D1550','A1550'])
    total_780  = sum(counts[key] for key in ['H780','V780','D780','A780'])
    total_CC   = sum(counts[key] for key in ['H1550H780','H1550V780','H1550D780','H1550A780','V1550H780','V1550V780','V1550D780','V1550A780','D1550H780','D1550V780','D1550D780','D1550A780','A1550H780','A1550V780','A1550D780','A1550A780'])
    accidentals = np.sum(np.prod(np.array(list(product([counts['H1550'],counts['V1550'],counts['D1550'],counts['A1550']], [counts['H780'],counts['V780'],counts['D780'],counts['A780']]))),axis=1))*Tcc

    intrinsic_heralding_1550 = -10*np.log10(total_CC/total_1550)-loss
    intrinsic_heralding_780 = -10*np.log10(total_CC/total_780)
    # The qber and qx will be overestimated at 30db
    qber = (counts['V1550H780']+counts['H1550V780'])/(counts['V1550H780']+counts['H1550V780']+counts['H1550H780']+counts['V1550V780'])
    qx = (counts['A1550D780']+counts['D1550A780'])/(counts['A1550D780']+counts['D1550A780']+counts['D1550D780']+counts['A1550A780'])

    Brightness = total_1550*total_780 / (total_CC - accidentals)
    Brightness /= power  # we want to normalise this for power.
    # Brightness *= 10**(-loss/10)  # we want to normalise this for loss.

    print(f'{Brightness=}, {qx=}, {qber=}')
    print(f'{intrinsic_heralding_1550=}, {intrinsic_heralding_780=}')
    f.close()
    return intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc

def plot_model(params, powers, loss_range,t_delta=0.5e-9,DC_A=100, DC_B=80, t_dead_A=25e-9, t_dead_B=45e-9):
    intrinsic_heralding_1550, intrinsic_heralding_780, qber, qx, Brightness, Tcc = params
    labels = ['intrinsic_heralding_1550: ', 'intrinsic_heralding_780: ','bit_err: ', 'phase_err: ', 'Brightness: ','t_delta: ']
    print([labels[i]+str(params[i]) for i in range(len(params))]+[f'DC_A={DC_A}, DC_B={DC_B}'])
    
    setup = secure_key_rates(d=4, t_delta=t_delta, DC_A=DC_A, DC_B=DC_B, e_b=qber, e_p=qx, t_dead_A=t_dead_A, t_dead_B=t_dead_B, loss_format='dB', custom=True)
    losses = loss_range
    for power in powers:
        akr = []
        for loss in losses:
            akr += [setup.custom_performance(tcc=Tcc, B=Brightness*power,eff_A=intrinsic_heralding_1550, eff_B=intrinsic_heralding_780+loss)]
           
        
        return akr

def plot_data(file_name, losses, power, time):
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
    # plot the overpass fitted from '/Users/thomasjaeken/Heriot-Watt University Team Dropbox/RES_EPS_EMQL/projects/Optical ground station/__experimental-work__/__BBM92__/loss profile for micius uplink/loss profile for micius uplink.csv'
    time = np.linspace(0,360,360)
    loss_behaviour = lambda x,freq,a,b,c,d: a*np.log10(b**2+c**2-2*b*c*np.cos(freq*(x-d))) - 8.9
    loss_overpass = partial(loss_behaviour,freq=-8.13786355e-05,  a=1.26519905e+01,  b=6.57313691e+03,  c=6.52804180e+03,d=1.82135169e+02)
    axloss = ax.twinx()
    axloss.plot(time,loss_overpass(time),'--', color='black')
    # axloss.ticklabel_format(axis='y', style='plain')
    # axloss.set_yscale('log')
    axloss.set_ylabel('Loss (dB)', color='black')
    # axloss.tick_params(axis='y', labelcolor='blue')

    # plot the 
    run_name = 'datarun_30_07_24'
    start_data, loss, power = find_start_point(run_name)
    params = [7.761328918344407, 7.542259886343475, 0.0335677812551474, 0.02531375770896907, 10826895.017621633, 4e-10]
    # params = [7.761328918344407, 7.542259886343475, 0.05435677812551474, 0.07731375770896907, 8826895.017621633, 4e-10]
    akr=plot_model(params, powers=[5], loss_range=loss_overpass(time),t_delta=0.5e-9,DC_A=200, DC_B=120, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(time,akr,color="darkblue")
    datapoints,yerr = plot_data(run_name,loss_overpass(time), 5, time)
    ax.errorbar(time[np.where(datapoints>0)],datapoints[np.where(datapoints>0)],yerr=yerr[np.where(datapoints>0)],linestyle='',ms=4,marker='o',color="darkblue", label="0km",markeredgecolor='black', markeredgewidth=1/2)
    # ax.scatter(time[np.where(datapoints>0)], datapoints[np.where(datapoints>0)], s=4,marker='o',  color="darkblue", edgecolors='black',linewidths=1/2)
    # ax.fill_between(time[np.where(datapoints!=0)], datapoints[np.where(datapoints!=0)]-yerr[np.where(datapoints!=0)], datapoints[np.where(datapoints!=0)]+yerr[np.where(datapoints!=0)],  color="darkblue",alpha=0.4)

    print(f"total key: {np.trapz(x=time[np.where(np.array(akr)>=0)],y=np.array(akr)[np.where(np.array(akr)>=0)]):.0f}bit, peak rate: {max(akr)}, nonzero time: {time[np.where(np.array(akr)>=0)][-1]-time[np.where(np.array(akr)>=0)][0]}" )

    # plot the 
    run_name = 'datarun_17_07_24'
    start_data, loss, power = find_start_point(run_name)
    #Brightness,t_delta, DC_A, DC_B, intrinsic_heralding_1550, intrinsic_heralding_780,bit_err, phase_err
    params = [9.47893858965945,9.89101987420694, 0.037228590019520497, 0.03834132110856429, 6135831.8248959305, 1e-09]
    #[9.89860936e-01, 3.42282186e-01, 2.00000000e+02, 6.99995871e+01,7.47866212e+00, 7.48347036e+00, 5.49766649e-02, 9.09056315e-02]
    akr=plot_model(params, powers=[5], loss_range=loss_overpass(time),t_delta=0.4e-9,DC_A=200, DC_B=70, t_dead_A=25e-9, t_dead_B=45e-9)
    ax.plot(time,akr,color="lightblue")
    datapoints,yerr = plot_data(run_name,loss_overpass(time), 5, time)
    ax.errorbar(time[np.where(datapoints>0)],datapoints[np.where(datapoints>0)],yerr=yerr[np.where(datapoints>0)],linestyle='',ms=4,marker='o',color="lightblue", label="10km",markeredgecolor='black', markeredgewidth=1/2)
    # ax.scatter(time[np.where(datapoints>0)], datapoints[np.where(datapoints>0)], s=4,marker='o',  color="lightblue", edgecolors='black',linewidths=1/2)
    # ax.fill_between(time[np.where(datapoints!=0)], datapoints[np.where(datapoints!=0)]-yerr[np.where(datapoints!=0)], datapoints[np.where(datapoints!=0)]+yerr[np.where(datapoints!=0)],  color="lightblue",alpha=0.4)
        
    print(f"total key: {np.trapz(x=time[np.where(np.array(akr)>=0)],y=np.array(akr)[np.where(np.array(akr)>=0)]):.0f}bit, peak rate: {max(akr)}, nonzero time: {time[np.where(np.array(akr)>=0)][-1]-time[np.where(np.array(akr)>=0)][0]}" )



    ax.legend(fontsize=6, loc=(0,0.4), frameon=False)#(0.8,.5))
    ax.set_xlabel('Time (s)',fontsize=10)
    ax.set_ylabel('Secret key rate (bits/s)',fontsize=10)#,color='blue')
    plt.tight_layout()
    plt.show()
