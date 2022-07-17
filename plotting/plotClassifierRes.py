import numpy as np
import matplotlib.pyplot as plt
import json

json_file = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/Classifier_stats.json'
save_fig_repo = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/plots/'


# oPs_scaled_logLdiff = oPs_logLdiff / oPs_nhits
# elec_scaled_logLdiff = elec_logLdiff / elec_nhits

############ Plotting Funcs #############

def plotnhits(nhits):
    plt.figure(figsize=(10, 7), dpi=100)
    bins = 10**(np.linspace(0, 4, 200))
    for particle in nhits:
        plt.hist(x = nhits[particle], bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('Nhits')
    plt.xscale('log')
    plt.yscale('log')
    # plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_nhits.png', dpi=1000)
    plt.close()

def plotLL(data):
    plt.figure(figsize=(10, 7), dpi=100)
    bins = np.arange(-0.02, 0.08, 0.001)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.hist(x = scaled_LLdiff, bins = bins, alpha=0.5, rwidth=0.85, label = particle)
    plt.legend(loc='best')
    #plt.xlim([-50, 50])
    plt.xlabel('log-Likelyhood difference divided by nhit')
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_LLdiff.png', dpi=1000)
    plt.close()

def plotLLvsDelay(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['delays'], scaled_LLdiff, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('decay time (ns)')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Decay Times')
    plt.legend(loc='best')
    #plt.xlim([-10, 100])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_LLdiff_delays.png', dpi=1000)
    plt.close()

def plotLLvsEnergy(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['energies'], scaled_LLdiff, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Energies')
    plt.legend(loc='best')
    plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_LLdiff_energies.png', dpi=1000)
    plt.close()

def plotLLvsNhit(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        scaled_LLdiff = data[particle]['classifications'] / data[particle]['nhits']
        plt.scatter(data[particle]['nhits'], scaled_LLdiff, label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('log-Likelyhood difference divided by nhit')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Energies')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_LLdiff_nhits.png', dpi=1000)
    plt.close()

def plotEnergyvsNhit(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['nhits'], data[particle]['energies'], label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('Reconstructed Event Energy (MeV)')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Energies')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_energies_nhits.png', dpi=1000)
    plt.close()

def plotDelayvsEnergy(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['energies'], data[particle]['delays'], label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Reconstructed Event Energy (MeV)')
    plt.ylabel('Delays (ns)')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Energies')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_delays_nhits.png', dpi=1000)
    plt.close()

def plotDelayvsNhit(data):
    plt.figure(figsize=(10, 7), dpi=100)
    for particle in data:
        plt.scatter(data[particle]['nhits'], data[particle]['delays'], label=particle)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Nhits')
    plt.ylabel('Delays (ns)')
    plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Energies')
    plt.legend(loc='best')
    #plt.xlim([-1, 5])
    #plt.ylim([-0.14, 0.065])
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_delays_nhits.png', dpi=1000)
    plt.close()

def plotROC(data):
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
    plt.title('Positronium Classifier ROC curve')
    #plt.show()
    plt.savefig(save_fig_repo + 'o-Ps_ROC_curve.png', dpi=1000)
    plt.close()

############### MAIN ################

def main():
    f = open(json_file)
    data = json.load(f)
    f.close()

    new_data = {}
    full_nhits = {}
    for particle in data:
        full_nhits[particle] = []
        new_data[particle] = {}
        new_data[particle]['classifications'] = []
        new_data[particle]['delays'] = []
        new_data[particle]['nhits'] = []
        new_data[particle]['energies'] = []
        for energy in data[particle]:
            if not (energy == '1.0' or (energy == '2.0' and particle == 'e-')):
                continue
            for i in range(len(data[particle][energy]['Classifier_results'])):
                recon_energy = data[particle][energy]['recon_event_energies'][i]
                nhits = data[particle][energy]['nhits'][i]
                full_nhits[particle].append(nhits)
                #if recon_energy > 50 or recon_energy < 1:
                if nhits > 1000 or nhits < 200:
                    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                    print(i)
                    print('particle: ', particle)
                    print('simulated particle energy: ', energy)
                    print('Classifier_result: ', data[particle][energy]['Classifier_results'][i])
                    print('delay: ', data[particle][energy]['delays'][i])
                    print('nhit: ', data[particle][energy]['nhits'][i])
                    print('recon_event_energy: ', data[particle][energy]['recon_event_energies'][i])
                    #continue
                new_data[particle]['classifications'].append(data[particle][energy]['Classifier_results'][i])
                new_data[particle]['delays'].append(data[particle][energy]['delays'][i])
                new_data[particle]['nhits'].append(nhits)
                new_data[particle]['energies'].append(recon_energy)
        new_data[particle]['classifications'] = np.array(new_data[particle]['classifications'])
        new_data[particle]['delays'] = np.array(new_data[particle]['delays'])
        new_data[particle]['nhits'] = np.array(new_data[particle]['nhits'])
        new_data[particle]['energies'] = np.array(new_data[particle]['energies'])

    # Plotting:
    plotnhits(full_nhits)
    plotLL(new_data)
    plotLLvsDelay(new_data)
    plotLLvsEnergy(new_data)
    plotLLvsNhit(new_data)
    plotEnergyvsNhit(new_data)
    plotDelayvsEnergy(new_data)
    plotDelayvsNhit(new_data)
    if 'o-Ps' in new_data and 'e-' in new_data:
        if len(new_data['o-Ps']['classifications']) > 0 and len(new_data['e-']['classifications']) > 0:
            plotROC(new_data)

if __name__ == '__main__':
    main()