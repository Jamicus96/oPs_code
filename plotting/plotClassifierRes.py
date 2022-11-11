import numpy as np
from numpy import unravel_index
import json

# json_file = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/Classifier_stats.json'
# save_fig_repo = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/plots/E_range/'
# show = False

# json_file = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/Classifier_stats.json'
# save_fig_repo = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/Plots/'
json_file = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/spectra/Classifier_stats.json'
save_fig_repo = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/spectra/plots/'
show = True

if not show:
    import matplotlib
    matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

######## OTHER #######

def applyCuts_old(data, energy_lim):
    new_data = {}
    full_nhits = {}
    for particle in data:
        if not (particle == 'o-Ps' or particle == 'e-'):
            continue
        full_nhits[particle] = []
        new_data[particle] = {}
        new_data[particle]['classifications'] = []
        new_data[particle]['delays'] = []
        new_data[particle]['nhits'] = []
        new_data[particle]['recon_event_energies'] = []
        for energy in data[particle]:
            if not (energy == '0.5-9.0' or energy == '1.5-10.0'):
                continue
            for i in range(len(data[particle][energy]['Classifier_results'])):
                recon_energy = data[particle][energy]['recon_event_energies'][i]
                nhits = data[particle][energy]['nhits'][i]
                full_nhits[particle].append(nhits)
                # Apply cut
                if nhits > 6000 or nhits < 200 or recon_energy < energy_lim[0] or recon_energy > energy_lim[1]:
                    continue
                new_data[particle]['classifications'].append(data[particle][energy]['Classifier_results'][i])
                new_data[particle]['delays'].append(data[particle][energy]['delays'][i])
                new_data[particle]['nhits'].append(nhits)
                new_data[particle]['recon_event_energies'].append(recon_energy)
        new_data[particle]['classifications'] = np.array(new_data[particle]['classifications'])
        new_data[particle]['delays'] = np.array(new_data[particle]['delays'])
        new_data[particle]['nhits'] = np.array(new_data[particle]['nhits'])
        new_data[particle]['recon_event_energies'] = np.array(new_data[particle]['recon_event_energies'])

    return new_data, full_nhits

cut_particles = ['IBD', 'alphaN13C']
cut_simEnergies = ['0.5-9.0']

def general_Cuts(data, new_data):
    for i in range(len(data['nhits'])):
        for key in data:
            new_data[key].append(data[key][i])
    return new_data

def IBD_Cuts(data, new_data):
    print('Doing IBD cuts...')

    for i in range(len(data['nhits'])):
        recon_energy = data['recon_event_energies'][i]
        nhits = data['nhits'][i]
        radius = data['radius'][i]
        alphaN_classifier = data['alphaNreactor_lassifier_results'][i]
        if nhits > 6000 or nhits < 130:
            continue
        if recon_energy > 3.5:
            continue
        if radius > 5700:  # in mm
            continue
        # if alphaN_classifier < 5:
        #     continue
        for key in data:
            new_data[key].append(data[key][i])
    return new_data

def alphaN13C_Cuts(data, new_data):
    '''Do coincidence tagging and other cuts'''
    print('Doing alphaN13C cuts...')
    # ● FV: 5.7 m
    # ● Prompt energy: 0.9 – 8.0 MeV
    # ● Delayed energy: 1.85 – 2.4 MeV
    # ● Delta R: 1.5 m
    # ● Delta t: 400 ns – 0.8 ms 
    previous_entryNum = -1
    for i in range(len(data['nhits'])):
        entry_num = data['entry_num'][i]
        if previous_entryNum != entry_num:
            prompt_energy = data['recon_event_energies'][i]
            prompt_nhits = data['nhits'][i]
            prompt_radius = data['radius'][i]
            prompt_time = data['delta_times'][i]
            prompt_alphaN_classifier = data['alphaNreactor_lassifier_results'][i]

            previous_entryNum = entry_num
        else:
            delayed_energy = data['recon_event_energies'][i]
            delayed_nhits = data['nhits'][i]
            delayed_radius = data['radius'][i]
            delayed_time = data['delta_times'][i]

            # Apply cuts
            #print('prompt_energy:', prompt_energy)
            if prompt_energy < 0.9 or prompt_energy > 8.0:
                continue
            #print('delayed_energy:', delayed_energy)
            if delayed_energy < 1.85 or delayed_energy > 2.4:
                continue
            #print('prompt_radius:', prompt_radius, 'delayed_radius:', delayed_radius)
            if prompt_radius > 5700 or delayed_radius > 5700:  # in mm
                continue
            if abs(delayed_radius - prompt_radius) > 1500:  # in mm  FIXME
                continue
            if delayed_time < 400 or delayed_time > 0.8E6:  # in ns
                continue
            if prompt_energy > 3.5:
                continue
            # if prompt_alphaN_classifier < 5:
            #     continue

            # Record data that passes cuts
            for key in data:
                # Record prompt event
                new_data[key].append(data[key][i - 1])
                # Do not record delayed event
                # new_data[key].append(data[key][i])

            previous_entryNum = entry_num
    return new_data

def applyCuts_new(data, energy_lim):
    cut_funcs = {
        'IBD': IBD_Cuts,
        'alphaN13C': alphaN13C_Cuts
    }

    new_data = {}
    full_nhits = {}
    for particle in data:
        if particle not in cut_particles:
            continue
        full_nhits[particle] = []
        new_data[particle] = {}
        for key in data[particle][cut_simEnergies[0]]:
            new_data[particle][key] = []

        for sim_energy in data[particle]:
            if sim_energy not in cut_simEnergies:
                continue
            full_nhits[particle].append(data[particle][sim_energy]['nhits'])
            # Apply cuts
            if particle in cut_funcs:
                new_data[particle] = cut_funcs[particle](data[particle][sim_energy], new_data[particle])
            else:
                new_data[particle] = general_Cuts(data[particle][sim_energy], new_data[particle])

        for key in new_data[particle]:
            new_data[particle][key] = np.array(new_data[particle][key])
            print('N_data[' + particle + '][' + cut_simEnergies[0] + '][' + key + '] =', len(data[particle][cut_simEnergies[0]][key]))
            print('N_data_cut[' + particle + '][' + cut_simEnergies[0] + '][' + key + '] =', len(new_data[particle][key]))

    return new_data, full_nhits

############ Plotting Funcs #############

classifier_name = 'alphaNreactor_lassifier_results'
# classifier_name = 'oPs_lassifier_results'

def plotLLvsEnergy(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['recon_event_energies'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('log-Likelihood difference')
    plt.title('Likelihood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_energies_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotnhits(nhits, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    bins = 10**(np.linspace(0, 4, 200))
    for particle in nhits:
        plt.hist(x = nhits[particle], bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('Nhits')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Nhits, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLL(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    bins = np.linspace(np.min(data['IBD'][classifier_name]), np.max(data['IBD'][classifier_name]), 100)
    N = len(data['IBD'][classifier_name])
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.hist(x = scaled_LLdiff[:N], bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('log-Likelihood difference')
    plt.title('LL difference, for reconstructed energy below 3.5 MeV')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotEnergy(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    bins = np.arange(1, 8.5, 0.1)
    N = len(data['IBD'][classifier_name])
    for particle in data:
        energies = data[particle]['recon_event_energies']# / data[particle]['nhits']
        print('Number of ' + particle + ' events: ' + str(len(energies)))
        plt.hist(x = energies[:N], bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Unitless')
    plt.title('Energy Spectrum')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_Energies_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsDelay(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['delays'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('decay time (ns)')
    plt.ylabel('log-Likelihood difference')
    plt.title('Likelihood o-Ps vs e+ For Different Decay Times, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsRadius(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['radius'] / 1000.0, scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('radius (m)')
    plt.ylabel('log-Likelihood difference')
    plt.title('Likelihood o-Ps vs e+ For Different Decay Times, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsdeltaT(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['delta_times'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('delta times (ns)')
    plt.ylabel('log-Likelihood difference')
    plt.title('Likelihood o-Ps vs e+ For Different Decay Times, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsEntryNum(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['entry_num'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('entry number')
    plt.ylabel('log-Likelihood difference')
    plt.title('Likelihood o-Ps vs e+ For Different Decay Times, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsNhit(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
        scaled_LLdiff = data[particle][classifier_name]
        plt.scatter(data[particle]['nhits'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('log-Likelihood difference divided by nhit')
    plt.title('Likelihood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotEnergyvsNhit(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['nhits'], data[particle]['recon_event_energies'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('Reconstructed Event Energy (MeV)')
    plt.title('Likelihood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_energies_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotDelayvsEnergy(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['recon_event_energies'], data[particle]['delays'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('Delays (ns)')
    plt.title('Likelihood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_delays_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotDelayvsNhit(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['nhits'], data[particle]['delays'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('Delays (ns)')
    plt.title('Likelihood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_delays_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsLL(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        nx = np.linspace(-40, 80, 100)
        ny = np.linspace(-20, 25, 100)
        # plt.hist2d(data[particle]['alphaNreactor_lassifier_results'], data[particle]['oPs_lassifier_results'], bins=(nx, ny), label = particle)
        plt.scatter(data[particle]['alphaNreactor_lassifier_results'], data[particle]['oPs_lassifier_results'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('alphaNreactor log-Likelihood difference')
    plt.ylabel('o-Ps log-Likelihood difference')
    plt.title('Likelihood Difference Conparison, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_energies_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotROC(data, energy_lim):
    # scaled_oPs_LLdiff = data['o-Ps']['classifications']# / data['o-Ps']['nhits']
    # scaled_elec_LLdiff = data['e-']['classifications']# / data['e-']['nhits']

    # scaled_LLdiff = data[particle]['classifications']# / data[particle]['nhits']
    scaled_IBD_LLdiff = data['IBD'][classifier_name]
    scaled_alphaN_LLdiff = data['alphaN13C'][classifier_name]

    N = 100
    logL_perNHit = np.linspace(-15, 20, N)
    signal = np.zeros(N)
    noise = np.zeros(N)

    for n in range(N):
        cut = logL_perNHit[n]

        # Ratio of signal (o-Ps) that passes the cut
        for i in range(len(scaled_IBD_LLdiff)):
            if scaled_IBD_LLdiff[i] > cut:
                signal[n] += 1
        signal[n] *= 100 / len(scaled_IBD_LLdiff)

        # Ratio of noise (e-) that passes the cut
        for j in range(len(scaled_alphaN_LLdiff)):
            if scaled_alphaN_LLdiff[j] > cut:
                noise[n] += 1
        noise[n] *= 100 / len(scaled_alphaN_LLdiff)

    plt.figure(figsize=(10, 7), dpi=100)
    plt.plot(noise, signal)
    plt.xlabel('Percentage of noise (reactor alphaN) that passes cut')
    plt.ylabel('Percentage of signal (IBD) that passes cut')
    plt.title('Positronium Classifier ROC curve, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_ROC_curve_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plot_mult_ROC_oPs(data, energy_lim):
    oPs_signal_LLdiff = data['IBD']['oPs_lassifier_results']
    oPs_noise_LLdiff = data['alphaN13C']['oPs_lassifier_results']
    alphaN_signal_LLdiff = data['IBD']['alphaNreactor_lassifier_results']
    alphaN_noise_LLdiff = data['alphaN13C']['alphaNreactor_lassifier_results']

    N = 50
    oPs_logL = np.linspace(-20, 30, N)
    alphaN_logL = np.array([-40, -20, -10, -5, 0, 5])

    for alphaN_cut in alphaN_logL:
        signal = np.zeros(N)
        noise = np.zeros(N)
        for m in range(N):
            oPs_cut = oPs_logL[m]

            # Ratio of signal (o-Ps) that passes the cut
            signal_tot = 0
            for i in range(len(oPs_signal_LLdiff)):
                if alphaN_signal_LLdiff[i] > alphaN_cut:
                    signal_tot +=1
                    if oPs_signal_LLdiff[i] > oPs_cut:
                        signal[m] += 1
            signal[m] *= 100 / signal_tot

            # Ratio of noise (e-) that passes the cut
            noise_tot = 0
            for j in range(len(oPs_noise_LLdiff)):
                if alphaN_noise_LLdiff[j] > alphaN_cut:
                    noise_tot += 1
                    if oPs_noise_LLdiff[j] > oPs_cut:
                        noise[m] += 1
            noise[m] *= 100 / noise_tot

        plt.plot(noise, signal, label = 'alphaN cut: ' + str(alphaN_cut))
    # plt.xlim((0, 100))
    # plt.ylim((0, 100))
    plt.legend(loc='best')
    plt.xlabel('noise (%)')
    plt.ylabel('signal (%)')
    plt.title('o-Ps ROC Curve for Different Reactor alpha-n Classifier Cuts')
    plt.show()

def plot_mult_ROC_alphaN(data, energy_lim):
    oPs_signal_LLdiff = data['IBD']['oPs_lassifier_results']
    oPs_noise_LLdiff = data['alphaN13C']['oPs_lassifier_results']
    alphaN_signal_LLdiff = data['IBD']['alphaNreactor_lassifier_results']
    alphaN_noise_LLdiff = data['alphaN13C']['alphaNreactor_lassifier_results']

    N = 50
    oPs_logL = np.array([-30, -20, -10, 0, 5])
    alphaN_logL = np.linspace(-50, 100, N)

    for oPs_cut in oPs_logL:
        signal = np.zeros(N)
        noise = np.zeros(N)
        for m in range(N):
            alphaN_cut = alphaN_logL[m]

            # Ratio of signal (o-Ps) that passes the cut
            signal_tot = 0
            for i in range(len(oPs_signal_LLdiff)):
                if oPs_signal_LLdiff[i] > oPs_cut:
                    signal_tot +=1
                    if alphaN_signal_LLdiff[i] > alphaN_cut:
                        signal[m] += 1
            signal[m] *= 100 / signal_tot

            # Ratio of noise (e-) that passes the cut
            noise_tot = 0
            for j in range(len(oPs_noise_LLdiff)):
                if oPs_noise_LLdiff[j] > oPs_cut:
                    noise_tot += 1
                    if alphaN_noise_LLdiff[j] > alphaN_cut:
                        noise[m] += 1
            noise[m] *= 100 / noise_tot

        plt.plot(noise, signal, label = 'o-Ps cut: ' + str(oPs_cut))
    # plt.xlim((0, 100))
    # plt.ylim((0, 100))
    plt.legend(loc='best')
    plt.xlabel('noise (%)')
    plt.ylabel('signal (%)')
    plt.title('alpha-n ROC Curve for Different Reactor o-Ps Classifier Cuts')
    plt.show()

def plot2dSig(data, energy_lim):
    oPs_signal_LLdiff = data['IBD']['oPs_lassifier_results']
    oPs_noise_LLdiff = data['alphaN13C']['oPs_lassifier_results']
    alphaN_signal_LLdiff = data['IBD']['alphaNreactor_lassifier_results']
    alphaN_noise_LLdiff = data['alphaN13C']['alphaNreactor_lassifier_results']

    N = 50
    oPs_logL = np.linspace(-20, 30, N)
    alphaN_logL = np.linspace(-40, 80, N)
    signal = np.zeros((N, N))
    noise = np.zeros((N, N))

    for n in range(N):
        oPs_cut = oPs_logL[n]
        for m in range(N):
            alphaN_cut = alphaN_logL[m]

            # Ratio of signal (o-Ps) that passes the cut
            for i in range(len(oPs_signal_LLdiff)):
                if oPs_signal_LLdiff[i] > oPs_cut and alphaN_signal_LLdiff[i] > alphaN_cut:
                    signal[n, m] += 1
            signal[n, m] *= 100 / len(oPs_signal_LLdiff)

            # Ratio of noise (e-) that passes the cut
            for j in range(len(oPs_noise_LLdiff)):
                if oPs_noise_LLdiff[j] > oPs_cut and alphaN_noise_LLdiff[j] > alphaN_cut:
                    noise[n, m] += 1
            noise[n, m] *= 100 / len(oPs_noise_LLdiff)

            if n == 0 and m == 15:
                print('n =', n, ', m=', m)
                print('oPs_cut =', oPs_cut, ', alphaN_cut =', alphaN_cut)
                print('z =', signal[n, m] / np.sqrt(signal[n, m] + noise[n, m]))
            if n == 15 and m == 0:
                print('n =', n, ', m=', m)
                print('oPs_cut =', oPs_cut, ', alphaN_cut =', alphaN_cut)
                print('z =', signal[n, m] / np.sqrt(signal[n, m] + noise[n, m]))

    print('plotting...')
    x, y = np.meshgrid(oPs_logL, alphaN_logL)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    z = signal / np.sqrt(signal + noise)
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    x_max = 0
    y_max = 0
    for i in range(len(z)):
        for j in range(len(z[i])):
            if not np.isnan(z[i, j]):
                if z[i, j] > z[x_max, y_max]:
                    x_max = i
                    y_max = j

    print(x_max)
    print(y_max)
    print(z[x_max, y_max])
    print('Max value at oPs cut =', oPs_logL[x_max], ', and alphaN cut =', alphaN_logL[y_max])

    plt.xlabel("o-Ps classifier cut")
    plt.ylabel("Reactor alpha-n classifier cut")
    ax.set_zlabel("signal / sqrt(signal + noise)")
    Title = 'Classifier Performance'
    plt.title(Title)
    plt.show()



############### MAIN ################

def main():
    f = open(json_file)
    data = json.load(f)
    f.close()

    # Upper energy limit
    # energy_lims = [(0.9, 2.5), (2.5, 4.0), (4.0, 5.5), (5.5, 7.0), (7.0, 8.5), (8.5, 10.0)]
    energy_lims = [(0.9, 10.0)]
    for energy_lim in energy_lims:
        cut_data, full_nhits = applyCuts_new(data, energy_lim)

        # Plotting:
        # plotnhits(full_nhits, energy_lim)
        # plotLLvsRadius(cut_data, energy_lim)
        # plotLLvsdeltaT(cut_data, energy_lim)
        # plotLLvsEntryNum(cut_data, energy_lim)
        # plotLLvsDelay(cut_data, energy_lim)
        # plotLLvsEnergy(cut_data, energy_lim)
        # plotLLvsNhit(cut_data, energy_lim)
        # plotEnergyvsNhit(cut_data, energy_lim)
        # plotDelayvsEnergy(cut_data, energy_lim)
        # plotDelayvsNhit(cut_data, energy_lim)
        # plotEnergy(cut_data, energy_lim)
        # plotLLvsLL(cut_data, energy_lim)
        plotLL(cut_data, energy_lim)
        # plotROC(cut_data, energy_lim)
        # plot_mult_ROC_oPs(cut_data, energy_lim)
        # plot_mult_ROC_alphaN(cut_data, energy_lim)
        # plot2dSig(cut_data, energy_lim)

if __name__ == '__main__':
    main()