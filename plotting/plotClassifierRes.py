import numpy as np
import json

json_file = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/Classifier_stats.json'
save_fig_repo = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/plots/E_range/'
show = False

# json_file = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/Classifier_stats.json'
# save_fig_repo = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/Plots/'
# show = True

if not show:
    import matplotlib
    matplotlib.use('pdf')
import matplotlib.pyplot as plt

######## OTHER #######

def applyCuts(data, energy_lim):
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
        new_data[particle]['energies'] = []
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
                new_data[particle]['energies'].append(recon_energy)
        new_data[particle]['classifications'] = np.array(new_data[particle]['classifications'])
        new_data[particle]['delays'] = np.array(new_data[particle]['delays'])
        new_data[particle]['nhits'] = np.array(new_data[particle]['nhits'])
        new_data[particle]['energies'] = np.array(new_data[particle]['energies'])

    return new_data, full_nhits

############ Plotting Funcs #############

def plotLLvsEnergy(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['energies'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
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
    bins = np.arange(-0.02, 0.03, 0.0007)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.hist(x = scaled_LLdiff, bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('log-Likelyhood difference divided by nhit')
    plt.title('LL difference, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsDelay(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['delays'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('decay time (ns)')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('Likelyhood o-Ps vs e+ For Different Decay Times, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsEnergy(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['energies'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_LLdiff_energies_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotLLvsNhit(data, energy_lim):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['nhits'], scaled_LLdiff, alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
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
        plt.scatter(data[particle]['nhits'], data[particle]['energies'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('Reconstructed Event Energy (MeV)')
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
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
        plt.scatter(data[particle]['energies'], data[particle]['delays'], alpha=0.5, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('Delays (ns)')
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
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
    plt.title('Likelyhood o-Ps vs e+ For Different Energies, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_delays_nhits_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

def plotROC(data, energy_lim):
    scaled_oPs_LLdiff = data['o-Ps']['classifications'] / data['o-Ps']['nhits']
    scaled_elec_LLdiff = data['e-']['classifications'] / data['e-']['nhits']

    N = 100
    logL_perNHit = np.linspace(-0.09, 0.05, N)
    signal = np.zeros(N)
    noise = np.zeros(N)

    for n in range(N):
        cut = logL_perNHit[n]

        # Ratio of signal (o-Ps) that passes the cut
        for i in range(len(scaled_oPs_LLdiff)):
            if scaled_oPs_LLdiff[i] > cut:
                signal[n] += 1
        signal[n] *= 100 / len(scaled_oPs_LLdiff)

        # Ratio of noise (e-) that passes the cut
        for j in range(len(scaled_elec_LLdiff)):
            if scaled_elec_LLdiff[j] > cut:
                noise[n] += 1
        noise[n] *= 100 / len(scaled_elec_LLdiff)

    plt.figure(figsize=(10, 7), dpi=100)
    plt.plot(noise, signal)
    plt.xlabel('Percentage of noise (e-) that passes cut')
    plt.ylabel('Percentage of signal (o-Ps) that passes cut')
    plt.title('Positronium Classifier ROC curve, for reconstructed energy in ' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + 'MeV')
    if show:
        plt.show()
    else:
        plt.savefig(save_fig_repo + 'o-Ps_ROC_curve_' + str(energy_lim[0]) + '-' + str(energy_lim[1]) + '.png', dpi=1000)
        plt.close()

############### MAIN ################

def main():
    f = open(json_file)
    data = json.load(f)
    f.close()

    # Upper energy limit
    #energy_lims = [(0.0, 2.5), (2.5, 4.0), (4.0, 5.5), (5.5, 7.0), (7.0, 8.5), (8.5, 10.0)]
    energy_lims = [(0.0, 12.0)]
    for energy_lim in energy_lims:
        cut_data, full_nhits = applyCuts(data, energy_lim)

        # Plotting:
        plotnhits(full_nhits, energy_lim)
        plotLLvsDelay(cut_data, energy_lim)
        plotLLvsEnergy(cut_data, energy_lim)
        plotLLvsNhit(cut_data, energy_lim)
        plotEnergyvsNhit(cut_data, energy_lim)
        plotDelayvsEnergy(cut_data, energy_lim)
        plotDelayvsNhit(cut_data, energy_lim)
        plotLL(cut_data, energy_lim)
        if 'o-Ps' in cut_data and 'e-' in cut_data:
            if len(cut_data['o-Ps']['classifications']) > 0 and len(cut_data['e-']['classifications']) > 0:
                plotROC(cut_data, energy_lim)

if __name__ == '__main__':
    main()