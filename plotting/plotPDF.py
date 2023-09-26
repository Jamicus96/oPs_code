import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import json

json_file = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/PDFs/PDF_stats.json'
save_fig_repo = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/PDFs/plots/'
datatype = 'Fitted'
#datatype = 'MC'
#datatype = 'both'
energies = True


f = open(json_file)
data = json.load(f)
f.close()

times = data['times']
elec_probs = data['e-']
oPs_probs = data['o-Ps']
# for energy in elec_probs:
#     plt.plot(times, elec_probs[energy]['PDF_' + datatype], ':', label='e-, ' + energy + 'MeV')
#     plt.plot(times, oPs_probs[energy]['PDF_' + datatype], '-', label='o-Ps, ' + energy + 'MeV')

if datatype == 'both':
    plt.figure(figsize=(10, 7), dpi=100)

    plt.plot(times, elec_probs['2.0']['PDF_MC'], ':', label='e-, 2.0MeV, MC')
    plt.plot(times, elec_probs['2.0']['PDF_Fitted'], ':', label='e-, 2.0MeV, Fitted')
    plt.plot(times, oPs_probs['1.0']['PDF_MC'], '-', label='o-Ps, 1.0MeV, MC')
    plt.plot(times, oPs_probs['1.0']['PDF_Fitted'], '-', label='o-Ps, 1.0MeV, Fitted')
    plt.plot(times, data['gamma']['2.0']['PDF_MC'], '--', label='gamma, 2.0MeV, MC')
    plt.plot(times, data['gamma']['2.0']['PDF_Fitted'], '--', label='gamma, 2.0MeV, Fitted')
    
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r'$t_{res}$ (ns)')
    plt.ylabel('probability density')
    plt.title('Pulse Shapes (~10 000 events each), MC and Fitted Positions')
    plt.xlim([-10, 200])
    plt.legend(loc='best')
    #plt.show()
    #plt.savefig(save_fig_repo + 'o-Ps_e-Energies' + datatype + '.png', dpi=1000)
    plt.savefig(save_fig_repo + 'particle_PDFs_both.png', dpi=1000)
    plt.close()

else:
    if energies:
        save_fig_repo += 'energies/'
        energy_list = [('0.5', '1.5', 'b'), ('3.0', '4.0', 'g'), ('5.0', '6.0', 'r'), ('7.0', '8.0', 'c'), ('9.0', '10.0', 'm')]
        plt.figure(figsize=(10, 7), dpi=100)
        for oPs_energy, elec_energy, colour in energy_list:
            plt.plot(times, elec_probs[elec_energy]['PDF_' + datatype], colour + ':', label='e-, ' + elec_energy + 'MeV')
            plt.plot(times, oPs_probs[oPs_energy]['PDF_' + datatype], colour + '-', label='o-Ps, ' + oPs_energy + 'MeV')

        #plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel(r'$t_{res}$ (ns)')
        plt.ylabel('probability density')
        plt.title('Pulse Shapes (~10 000 events each), ' + datatype + ' Position')
        plt.xlim([-10, 200])
        plt.legend(loc='best')
        #plt.show()
        #plt.savefig(save_fig_repo + 'o-Ps_e-Energies' + datatype + '.png', dpi=1000)
        plt.savefig(save_fig_repo + 'particle_PDFs_energies_' + datatype + '.png', dpi=1000)
        plt.close()

    else:
        plt.figure(figsize=(10, 7), dpi=100)

        plt.plot(times, elec_probs['2.0']['PDF_' + datatype], ':', label='e-, 2.0MeV')
        plt.plot(times, oPs_probs['1.0']['PDF_' + datatype], '-', label='o-Ps, 1.0MeV')
        plt.plot(times, data['gamma']['2.0']['PDF_' + datatype], '-', label='gamma, 2.0MeV')

        #plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel(r'$t_{res}$ (ns)')
        plt.ylabel('probability density')
        plt.title('Pulse Shapes (~10 000 events each), ' + datatype + ' Position')
        plt.xlim([-10, 200])
        plt.legend(loc='best')
        #plt.show()
        #plt.savefig(save_fig_repo + 'o-Ps_e-Energies' + datatype + '.png', dpi=1000)
        plt.savefig(save_fig_repo + 'particle_PDFs_' + datatype + '.png', dpi=1000)
        plt.close()