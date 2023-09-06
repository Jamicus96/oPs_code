import argparse
import numpy as np
import matplotlib.pyplot as plt

def read_info(in_file):
    # Open file
    print('Reading in:', in_file)
    f = open(in_file)
    lines = f.readlines()
    f.close()

    # Create time bin centers and time bin limits from first line
    time_lims = lines[0].split(', ')  # Nbins, lower_lim, upper_lim
    Nbins = int(time_lims[0])
    lower_lim = float(time_lims[1])
    upper_lim = float(time_lims[2])
    step = (upper_lim - lower_lim) / float(Nbins)
    time_bin_centers = np.linspace(lower_lim + (0.5 * step), upper_lim - (0.5 * step), Nbins)
    time_bin_lims = np.linspace(lower_lim, upper_lim, Nbins + 1)

    # Loop through the rest of the lines to construct the data tree.
    info = {}
    info = {}
    info['entry'] = []; info['evt'] = []; info['num_evts'] = []; info['max_evt'] = []
    info['Class'] = []; info['E'] = []; info['pos'] = []
    info['Delta_T'] = []; info['Nhit'] = []; info['t_res'] = []
    
    last_MC_entry, last_recon_entry = -1, -1
    num_evts_MC, num_evts_recon = 0, 0
    for i in range(1, len(lines)):
        line = lines[i].replace('\n', '').split(' ')

        if False:
            print('line =', line)

        if line[0] == 'MC':
            continue
        
        elif line[0] == 'recon':
            info['entry'].append(int(line[1])); info['evt'].append(int(line[2]))
            info['Class'].append(float(line[3])); info['E'].append(float(line[4]))
            info['pos'].append(line[5].split(',')); info['Delta_T'].append(float(line[6]))
            info['Nhit'].append(float(line[7]))
            info['t_res'].append(line[8].split(','))

            if int(line[1]) == last_recon_entry:
                info['num_evts'].append(num_evts_recon)
                info['max_evt'].append(int(line[2]))
                num_evts_recon += 1
                for j in range(1, num_evts_recon + 1):
                    info['num_evts'][-j] = num_evts_recon
                    info['max_evt'][-j] = int(line[2])
            else:
                last_recon_entry = int(line[1])
                num_evts_recon = 1
                info['num_evts'].append(num_evts_recon)
                info['max_evt'].append(int(line[2]))

        else:
            print('[WARNING] Line {} has unknown first element: {}'.format(i, line[0]))

    # Convert lists to numpy arrays
    for key_1 in info:
        if isinstance(info[key_1], list):
            if isinstance(info[key_1][0], list):
                info[key_1] = np.asarray(info[key_1], dtype=float)
            else:
                info[key_1] = np.array(info[key_1])
        else:
            print('[WARNING] info[{}][{}] is type {}, not list.'.format(key_1, type(info[key_1])))

    return info, time_bin_centers, time_bin_lims

def create_prompt_delay(info):

    # Create new data structure: keep only prompt events)
    re_info = {}

    # Find prompt and delayed event indices
    indices = np.intersect1d(np.intersect1d(np.where(info['evt'] == 0)[0],\
                                            np.where(info['num_evts'] == 2)[0], assume_unique=True),\
                                            np.where(info['max_evt'] == 1)[0], assume_unique=True)

    # Populate new data structure based on this
    for key_1 in info:
        if key_1 == 'evt' or key_1 == 'num_evts' or key_1 == 'max_evt':
            continue
        re_info[key_1] = info[key_1][indices]

    return re_info


def fit_line(x, y):
    '''Fit data to straight line with y(0) = 0'''
    return np.sum(y * x) / np.sum(x**2)

# Read info from files
IBD_file_info, IBD_times, IBD_time_bins = read_info('/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_6p0_TeDiol_1p5_bismsb_dda_scintillator_directors_review/getPDFs/info_ReactorIBD_0.txt')
alphaN_file_info, alphaN_times, alphaN_time_bins = read_info('/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_6p0_TeDiol_1p5_bismsb_dda_scintillator_directors_review/getPDFs/info_alphaN_13C_0.txt')

# Keep only prompt events
IBD_file_info = create_prompt_delay(IBD_file_info)
alphaN_file_info = create_prompt_delay(alphaN_file_info)

# Only keep events in energy range
E_max = 7.0
E_min = 2.0
# E_max = 9.5
# E_min = 0.5
idx_IBD = np.intersect1d(np.where(IBD_file_info['E'] <= E_max)[0],
                         np.where(IBD_file_info['E'] >= E_min)[0], assume_unique=True)
idx_alphaN = np.intersect1d(np.where(alphaN_file_info['E'] <= E_max)[0],
                            np.where(alphaN_file_info['E'] >= E_min)[0], assume_unique=True)
E_IBD = IBD_file_info['E'][idx_IBD]
E_alphaN = alphaN_file_info['E'][idx_alphaN]
Nhit_IBD = IBD_file_info['Nhit'][idx_IBD]
Nhit_alphaN = alphaN_file_info['Nhit'][idx_alphaN]

# Compute Nhits vs E slopes
IBD_slope = fit_line(E_IBD, Nhit_IBD)
print('IBD_slope =', IBD_slope)
alphaN_slope = fit_line(E_alphaN, Nhit_alphaN)
print('alphaN_slope =', alphaN_slope)
E = np.linspace(E_min, E_max, 1000)

# Plotting!
fig, (ax1, ax2) = plt.subplots(2, sharex=True, height_ratios=[3, 1])
fig.suptitle(r'Comparing Recon Energy to Nhits, for Reactor IDB and $\alpha$-n Events')

ax1.scatter(E_IBD, Nhit_IBD, marker='+', label='IBD')
ax1.scatter(E_alphaN, Nhit_alphaN, marker='+', label=r'$\alpha$-n')
ax1.vlines([3.5, 5.4], np.min(np.concatenate((Nhit_alphaN, Nhit_IBD), axis=0)),
           [0.5 * (IBD_slope + alphaN_slope) * 3.5, 0.5 * (IBD_slope + alphaN_slope) * 5.4],
           colors='k', linestyles='dashed')
ax1.hlines([0.5 * (IBD_slope + alphaN_slope) * 3.5, 0.5 * (IBD_slope + alphaN_slope) * 5.4],
           E_min, [3.5, 5.4], colors='k', linestyles='dashed')
print('Nhits 1 =', 0.5 * (IBD_slope + alphaN_slope) * 3.5)
print('Nhits 2 =', 0.5 * (IBD_slope + alphaN_slope) * 5.4)

ax1.plot(E, IBD_slope * E, label='IBD linear fit', color='b')
ax1.plot(E, alphaN_slope * E, label=r'$\alpha$-n linear fit', color='r')

ax1.set_ylabel('Nhits (cleaned)')
ax1.set_xlim((E_min, E_max))
ax1.set_ylim((np.min(np.concatenate((Nhit_alphaN, Nhit_IBD), axis=0)), np.max(np.concatenate((Nhit_alphaN, Nhit_IBD), axis=0))))
ax1.legend(loc='best')

counts_IBD, bins_IBD, bars_IBD = ax2.hist(E_IBD, label='IBD', bins=60)
counts_alphaN, bins_alphaN, bars_alphaN = ax2.hist(E_alphaN, label=r'$\alpha$-n', bins=60)
ax2.vlines([3.5, 5.4], np.min(np.concatenate((counts_IBD, counts_alphaN), axis=0)),
                       np.max(np.concatenate((counts_IBD, counts_alphaN), axis=0)),
                       colors='k', linestyles='dashed')
ax2.set_xlabel('Energy (MeV)')
ax2.set_xlim((E_min, E_max))

plt.show()
