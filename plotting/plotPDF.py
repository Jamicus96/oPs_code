import numpy as np
import matplotlib.pyplot as plt
import json

json_file = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/PDFs/PDF_stats.json'
save_fig_repo = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/PDFs/plots/'
#datatype = 'Fitted'
datatype = 'MC'


f = open(json_file)
data = json.load(f)
f.close()

times = data['times']
elec_probs = data['e-']
oPs_probs = data['o-Ps']
plt.figure(figsize=(10, 7), dpi=100)
for energy in elec_probs:
    plt.plot(times, elec_probs[energy]['PDF_' + datatype], ':', label='e-, ' + energy + 'MeV')
    plt.plot(times, oPs_probs[energy]['PDF_' + datatype], '-', label='o-Ps, ' + energy + 'MeV')

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$t_{res}$ (ns)')
plt.ylabel('probability density')
plt.title('Pulse Shapes (~10 000 events each), ' + datatype + ' Position')
plt.xlim([-10, 200])
plt.legend(loc='best')
#plt.show()
plt.savefig(save_fig_repo + 'o-Ps_e-Energies' + datatype + '.png', dpi=1000)
plt.close()