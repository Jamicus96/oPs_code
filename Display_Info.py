import tkinter as tk
import ttk
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

####################  READING FILE FUNCTIONS  ####################

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--verbose', '-v', type=bool, dest='verbose',
                        default=False, help='Print extra statements')
    parser.add_argument('-info_files', '-f', type=list, dest='info_files', help='List of text files with evtn info, in following order:\n\
                        [IBD_oldPDF_wo_oPs, alphaN_oldPDF_wo_oPs, IBD_oldPDF_wi_oPs, alphaN_oldPDF_wi_oPs,\
                        IBD_newPDF_wo_oPs, alphaN_newPDF_wo_oPs, IBD_newPDF_wi_oPs, alphaN_newPDF_wi_oPs]',
                        # default=[
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_IBD/oldPDF_without_oPs/info_Positronium_partial_IBD_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_alphaN/oldPDF_without_oPs/info_Positronium_partial_alphaN13C_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_IBD/oldPDF_with_oPs/info_Positronium_partial_IBD_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_alphaN/oldPDF_with_oPs/info_Positronium_partial_alphaN13C_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_IBD/newPDF_without_oPs/info_Positronium_partial_IBD_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_alphaN/newPDF_without_oPs/info_Positronium_partial_alphaN13C_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_IBD/newPDF_with_oPs/info_Positronium_partial_IBD_0.5-9.0.txt',
                        # '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/labppo_0p5_scintillator/PDFs_alphaN/newPDF_with_oPs/info_Positronium_partial_alphaN13C_0.5-9.0.txt',
                        # ])
                        default=[
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/info_Positronium_IBD_0.5-9.0_0.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/info_Positronium_alphaN13C_0.5-9.0_0.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/2p2labppo/nada.txt',
                        ])
    parser.add_argument('--ratdb_file', '-rf', type=str, dest='ratdb_file', help='ratdb file with old PDFs.',
                        default='/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/'\
                              + '2p2labppo/CLASSIFIER_ALPHAN_IBD_LIKELIHOOD.ratdb')
    parser.add_argument('--index', '-i', type=str, dest='index', help='Index to get PDFs from in ratdb file',
                        default='labppo_2p2_scintillator')
    parser.add_argument('--new_file', '-n', type=str, dest='new_file', help='New ratdb file to print PDFs to.',
                        default='/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/'\
                              + '2p2labppo/CLASSIFIER_ALPHAN_IBD_LIKELIHOOD_(new).ratdb')

    args = parser.parse_args()
    return args


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
    # info['MC'] = {}
    # info['MC']['entry'] = []; info['MC']['evt'] = []; info['MC']['num_evts'] = []; info['MC']['max_evt'] = []
    # info['MC']['Class'] = []; info['MC']['KE'] = []; info['MC']['PDG'] = []
    # info['MC']['pos'] = []; info['MC']['Delta_T'] = []; info['MC']['t_res'] = []
    info['recon'] = {}
    info['recon']['entry'] = []; info['recon']['evt'] = []; info['recon']['num_evts'] = []; info['recon']['max_evt'] = []
    info['recon']['Class'] = []; info['recon']['E'] = []; info['recon']['pos'] = []
    info['recon']['Delta_T'] = []; info['recon']['t_res'] = []
    
    last_MC_entry, last_recon_entry = -1, -1
    num_evts_MC, num_evts_recon = 0, 0
    for i in range(1, len(lines)):
        line = lines[i].replace('\n', '').split(' ')

        if False:
            print('line =', line)

        if line[0] == 'MC':
            continue
            # info['MC']['entry'].append(int(line[1])); info['MC']['evt'].append(int(line[2]))
            # # info['MC']['Class'].append(float(line[3])); info['MC']['KE'].append(float(line[4]))
            # info['MC']['Class'].append(0.0); info['MC']['KE'].append(float(line[3]))
            # info['MC']['PDG'].append(int(float(line[4]))); info['MC']['pos'].append(line[5].split(','));
            # info['MC']['Delta_T'].append(float(line[6])); info['MC']['t_res'].append(line[7].split(','))

            # if int(line[1]) == last_MC_entry:
            #     info['MC']['num_evts'].append(num_evts_MC)
            #     info['MC']['max_evt'].append(int(line[2]))
            #     num_evts_MC += 1
            #     for j in range(1, num_evts_MC + 1):
            #         info['MC']['num_evts'][-j] = num_evts_MC
            #         info['MC']['max_evt'][-j] = int(line[2])
            # else:
            #     last_MC_entry = int(line[1])
            #     num_evts_MC = 1
            #     info['MC']['num_evts'].append(num_evts_MC)
            #     info['MC']['max_evt'].append(int(line[2]))
        
        elif line[0] == 'recon':
            info['recon']['entry'].append(int(line[1])); info['recon']['evt'].append(int(line[2]))
            # info['recon']['Class'].append(float(line[3])); info['recon']['E'].append(float(line[4]))
            info['recon']['Class'].append(0.0); info['recon']['E'].append(float(line[3]))
            info['recon']['pos'].append(line[4].split(',')); info['recon']['Delta_T'].append(float(line[5]))
            info['recon']['t_res'].append(line[6].split(','))

            if int(line[1]) == last_recon_entry:
                info['recon']['num_evts'].append(num_evts_recon)
                info['recon']['max_evt'].append(int(line[2]))
                num_evts_recon += 1
                for j in range(1, num_evts_recon + 1):
                    info['recon']['num_evts'][-j] = num_evts_recon
                    info['recon']['max_evt'][-j] = int(line[2])
            else:
                last_recon_entry = int(line[1])
                num_evts_recon = 1
                info['recon']['num_evts'].append(num_evts_recon)
                info['recon']['max_evt'].append(int(line[2]))

        else:
            print('[WARNING] Line {} has unknown first element: {}'.format(i, line[0]))

    # Convert lists to numpy arrays
    for key_1 in info:
        for key_2 in info[key_1]:
            if isinstance(info[key_1][key_2], list):
                if isinstance(info[key_1][key_2][0], list):
                    info[key_1][key_2] = np.asarray(info[key_1][key_2], dtype=float)
                else:
                    info[key_1][key_2] = np.array(info[key_1][key_2])
            else:
                print('[WARNING] info[{}][{}] is type {}, not list.'.format(key_1, key_2, type(info[key_1][key_2])))

    return info, time_bin_centers, time_bin_lims

def create_prompt_delay(info):

    # Create new data structure: separating into prompt and delayed events (rest gets thrown out)
    re_info = {}
    # re_info['MC'] = {}; re_info['MC']['prompt'] = {}; re_info['MC']['delayed'] = {}
    re_info['recon'] = {}; re_info['recon']['prompt'] = {}; re_info['recon']['delayed'] = {}

    # Find prompt and delayed event indices
    indices = {}
    # indices['MC'] = {}
    indices['recon'] = {}
    # indices['MC']['prompt'] = np.intersect1d(np.intersect1d(np.where(info['MC']['evt'] == 0)[0],\
    #                                                         np.where(info['MC']['num_evts'] == 2)[0], assume_unique=True),\
    #                                                         np.where(info['MC']['max_evt'] == 1)[0], assume_unique=True)
    # indices['MC']['delayed'] = np.intersect1d(np.intersect1d(np.where(info['MC']['evt'] == 1)[0],\
    #                                                          np.where(info['MC']['num_evts'] == 2)[0], assume_unique=True),\
    #                                                          np.where(info['MC']['max_evt'] == 1)[0], assume_unique=True)

    indices['recon']['prompt'] = np.intersect1d(np.intersect1d(np.where(info['recon']['evt'] == 0)[0],\
                                                               np.where(info['recon']['num_evts'] == 2)[0], assume_unique=True),\
                                                               np.where(info['recon']['max_evt'] == 1)[0], assume_unique=True)
    indices['recon']['delayed'] = np.intersect1d(np.intersect1d(np.where(info['recon']['evt'] == 1)[0],\
                                                                np.where(info['recon']['num_evts'] == 2)[0], assume_unique=True),\
                                                                np.where(info['recon']['max_evt'] == 1)[0], assume_unique=True)

    # Populate new data structure based on this
    for key_1 in info:
        for key_2 in info[key_1]:
            if key_2 == 'evt' or key_2 == 'num_evts' or key_2 == 'max_evt':
                continue
            re_info[key_1]['prompt'][key_2] = info[key_1][key_2][indices[key_1]['prompt']]
            re_info[key_1]['delayed'][key_2] = info[key_1][key_2][indices[key_1]['delayed']]

    return re_info

def read_ratdb_file(file_address, index):
    # Open file
    f = open(file_address)
    lines = f.readlines()
    f.close()

    in_correct_table = False
    old_PDFs = {}
    old_PDFs['IBD'] = [[], [], []]
    old_PDFs['alpha-n'] = [[], [], []]
    for process in old_PDFs:
        if process == 'IBD':
            process_str = 'ibd'
        else:
            process_str = 'alpha'
        for line in lines:
            if 'index:' in line:
                if index in line:
                    in_correct_table = True
                else:
                    in_correct_table = False
            elif in_correct_table:
                if 'times:' in line:
                    times = line.replace('times:', '').replace(' ', '').replace('[', '').replace(']', '').replace('\n', '').split(',')
                    if times[-1] == '':
                        times = times[:-1]
                    times = np.array(times, dtype=float)
                else:
                    for i in range(3):
                        string = process_str + '_probability_' + str(i) + ':'
                        if string in line:
                            old_PDFs[process][i] = line.replace(string, '').replace(' ', '').replace('[', '').replace(']', '').replace('\n', '').split(',')
                            if old_PDFs[process][i][-1] == '':
                                old_PDFs[process][i] = old_PDFs[process][i][:-1]
                            old_PDFs[process][i] = np.array(old_PDFs[process][i], dtype=float)

    bin_width = times[1] - times[0]
    time_bins = np.append(times - (0.5 * bin_width), times[-1] + (0.5 * bin_width))

    return time_bins, times, old_PDFs


####################  CUT FUNCTIONS  ####################

def get_cuts(cuts_entries):
    cuts = {}
    for key in cuts_entries:
        cuts[key] = []
        for entry in cuts_entries[key]:
            cut = entry.get().replace(' ', '')
            if cut == '':
                cuts[key].append(False)
            else:
                if key == 'Both R' or key == 'Delta R' or key == 'Both Z':
                    # Convert m (input) to mm (saved evt values)
                    cuts[key].append(float(cut) * 1000.0)
                else:
                    cuts[key].append(float(cut))
        print('cuts[{}] = {}'.format(key, cuts[key]))
    return cuts

def apply_cut(values, cuts):
    '''Need only one set of values (for example energies E), and one set
    of cuts (for example cuts['Prompt E'] = [min, max]'''

    if cuts[0]:
        cut_min_idx = np.where(values > cuts[0])[0]
    else:
        cut_min_idx = np.arange(0, values.size, 1)
    # print('cut_min_idx.size =', cut_min_idx.size)

    if cuts[1]:
        cut_max_idx = np.where(values < cuts[1])[0]
    else:
        cut_max_idx = np.arange(0, values.size, 1)
    # print('cut_max_idx.size =', cut_max_idx.size)

    cut_idx = np.intersect1d(cut_min_idx, cut_max_idx, assume_unique=True)
    # print('cut_idx.size =', cut_idx.size)

    return cut_idx

def apply_global_cuts(cuts_entries, input_data):
    '''Should pass either recon, or MC dataset, not both. (both should have same structue, i.e. E not KE).
    Should also pass cuts_entries['global'] rather than just cuts_entries.'''

    # print('Init_len =', input_data['prompt']['E'].size)

    # Get cuts from cut entries
    cuts = get_cuts(cuts_entries)

    # Prepare values
    prompt_R = np.sqrt(np.sum(input_data['prompt']['pos']**2, axis=1))
    delayed_R = np.sqrt(np.sum(input_data['delayed']['pos']**2, axis=1))
    Delta_R = np.sqrt(np.sum((input_data['delayed']['pos'] - input_data['prompt']['pos'])**2, axis=1))
    

    # Get indices that pass each cut
    prompt_E_idx = apply_cut(input_data['prompt']['E'], cuts['Prompt E'])
    delayed_E_idx = apply_cut(input_data['delayed']['E'], cuts['Delayed E'])
    prompt_R_idx = apply_cut(prompt_R, cuts['Both R'])
    delayed_R_idx = apply_cut(delayed_R, cuts['Both R'])
    prompt_Z_idx = apply_cut(input_data['prompt']['pos'][:,2], cuts['Both Z'])
    delayed_Z_idx = apply_cut(input_data['delayed']['pos'][:,2], cuts['Both Z'])
    Delta_R_idx = apply_cut(Delta_R, cuts['Delta R'])
    Delta_T_idx = apply_cut(input_data['delayed']['Delta_T'], cuts['Delta T'])
    
    # Combine cuts
    cut_results = [delayed_E_idx, prompt_R_idx, delayed_R_idx, prompt_Z_idx, delayed_Z_idx, Delta_R_idx, Delta_T_idx]
    cut_idx = prompt_E_idx
    # print('Combine: cut_idx.size =', cut_idx.size)
    for cut_result in cut_results:
        cut_idx = np.intersect1d(cut_idx, cut_result, assume_unique=True)
        # print('Combine: cut_idx.size =', cut_idx.size)

    # Create new cut dataset
    output_data = {}
    for key_1 in input_data:
        output_data[key_1] = {}
        for key_2 in input_data[key_1]:
            output_data[key_1][key_2] = input_data[key_1][key_2][cut_idx]
            output_data[key_1][key_2] = input_data[key_1][key_2][cut_idx]

    return output_data

def apply_PDF_specific_cuts(cuts_entries, input_data):
    '''Should pass either recon, or MC dataset, not both. (both should have same structue, i.e. E not KE).
    Should also pass cuts_entries['specific'] rather than just cuts_entries.'''

    # Get cuts from cut entries
    cuts = get_cuts(cuts_entries)  

    # Get indices that pass each cut
    PDF_1_idx = apply_cut(input_data['prompt']['E'], cuts['PDF 1'])
    PDF_2_idx = apply_cut(input_data['prompt']['E'], cuts['PDF 2'])
    PDF_3_idx = apply_cut(input_data['prompt']['E'], cuts['PDF 3'])

    # Apply cuts to info
    PDF_cut_indices = [PDF_1_idx, PDF_2_idx, PDF_3_idx]

    # Create new cut dataset
    output_data = []
    for PDF_cut_idx in PDF_cut_indices:
        output_data_temp = {}
        for key_1 in input_data:
            output_data_temp[key_1] = {}
            for key_2 in input_data[key_1]:
                output_data_temp[key_1][key_2] = input_data[key_1][key_2][PDF_cut_idx]
                output_data_temp[key_1][key_2] = input_data[key_1][key_2][PDF_cut_idx]
        output_data.append(output_data_temp)


    return output_data


####################  PLOTTING FUNCTIONS  ####################

def plotting(fig, canvas, plotting_info, inputs, info, old_info):
    '''Choose which thing to plot, and call appropriate function.'''

    # Unpack information
    cuts_entries, checkboxes = inputs

    if checkboxes['which_plot'].get() == 'PDF':
        plot_t_res(fig, canvas, plotting_info, inputs, info, old_info)
    elif checkboxes['which_plot'].get() == 'Class':
        plot_alphaN_classifier(fig, canvas, plotting_info, inputs, info)
    else:
        print("ERROR: checkboxes['which_plot'] should be 'PDF' or 'Class', not '{}'".format(checkboxes['which_plot']))

def plot_t_res(fig, canvas, plotting_info, inputs, info, old_info):
    global PDFs

    # Unpack information
    cuts_entries, checkboxes = inputs
    entry_x_min, entry_x_max, Nbins, var_y_log = plotting_info
    old_time_bins, old_times, old_PDFs = old_info

    ################ Setup ################

    # Clear figure
    fig.clear()

    # adding the subplot
    plot1 = fig.add_subplot(111)

    # Set labels
    plot1.set_xlabel(r'$t_{res}$ (ns)')
    plot1.set_title('Normalised Time Residual')

    # Data type
    info_type = 'recon'

    colours = {
        'IBD' : ['r', 'tab:orange', 'y', 'tab:purple', 'tab:pink'],
        'alpha-n': ['b', 'tab:blue', 'tab:cyan', 'tab:gray', 'k']
    }
    linestyles = {
        'without_oPs': ':',
        'with_oPs': '-'
    }
    x_min = old_time_bins[0]
    x_max = old_time_bins[-1]


    ################ plotting ################

    PDFs = {}

    for process in checkboxes['process']:
        if not checkboxes['process'][process].get():
            continue
        for plot_name in checkboxes['PDF']:
            if not checkboxes['PDF'][plot_name].get():
                continue

            ################ plot old ratdb PDFs ################
            if plot_name == 'ratdb':
                for i, old_PDF in enumerate(old_PDFs[process]):
                    if checkboxes['cuts']['PDF' + str(i+1)].get():
                        # Plot old (ratdb) PDF
                        plot1.stairs(old_PDF, old_time_bins, linestyle='--', color=colours[process][i],\
                                     label = process + ', old (ratdb) PDF {}'.format(i+1))

            ################ plotting new PDFs ################
            else:
                label = process + ' ' + plot_name
                linestyle = linestyles[plot_name]
                if plot_name == 'without_oPs':
                    plot_name = process + '_oldPDF_wo_oPs'
                else:
                    plot_name = process + '_oldPDF_wi_oPs'

                # Upack info
                file_info, times, time_bin_lims = info[plot_name]
                file_info = file_info[info_type]
                bin_width = time_bin_lims[1] - time_bin_lims[0]

                if time_bin_lims[0] < x_min:
                    x_min = time_bin_lims[0]
                if time_bin_lims[-1] > x_max:
                    x_max = time_bin_lims[-1]
            
                ################ plotting the graph with no cuts ################

                if checkboxes['cuts']['noCUTS'].get():
                    # Sum all event time residuals, and remove underflow and overflow bins
                    tot_t_res = np.sum(file_info['prompt']['t_res'], axis=0)[1:-1]
                    # Normalise
                    tot_t_res_PDF = tot_t_res / (np.sum(tot_t_res) * bin_width)

                    # Plot
                    plot1.stairs(tot_t_res_PDF, time_bin_lims, color=colours[process][3], linestyle = linestyle,
                                 label = label + ', no cuts (N_evts = ' + str(file_info['prompt']['E'].size) + str(')'))

                ################ plot with general cuts ################
                
                glob_cut_info = apply_global_cuts(cuts_entries['global'], file_info)

                if checkboxes['cuts']['glbCUTS'].get():
                    # Sum all event time residuals
                    glob_cut_tot_t_res = np.sum(glob_cut_info['prompt']['t_res'], axis=0)[1:-1]
                    # Normalise
                    glob_cut_tot_t_res_PDF = glob_cut_tot_t_res / (np.sum(glob_cut_tot_t_res) * bin_width)

                    # Plot
                    plot1.stairs(glob_cut_tot_t_res_PDF, time_bin_lims, color=colours[process][4], linestyle = linestyle,
                                 label = label + ', global cuts (N_evts = ' + str(glob_cut_info['prompt']['E'].size) + str(')'))

                ################ plot with PDF-specific cuts ################

                spec_cut_info = apply_PDF_specific_cuts(cuts_entries['specific'], glob_cut_info)

                PDFs_temp = []
                for i, PDF_info in enumerate(spec_cut_info):
                    # Sum all event time residuals, then normalise
                    PDF_cut_tot_t_res = np.sum(PDF_info['prompt']['t_res'], axis=0)[1:-1]
                    # Normalise
                    PDFs_temp.append(PDF_cut_tot_t_res / (np.sum(PDF_cut_tot_t_res) * bin_width))
                    if checkboxes['cuts']['PDF' + str(i+1)].get():
                        # Plot
                        plot1.stairs(PDFs_temp[i], time_bin_lims, linestyle = linestyle, color=colours[process][i],\
                                    label = label + ', PDF {} (N_evts = {})'.format(i+1, PDF_info['prompt']['E'].size))

                if PDFs_temp != []:
                    if process == 'alpha-n':
                        PDFs['alpha_probability_'] = (plot_name, PDFs_temp)
                    else:
                        PDFs['ibd_probability_'] = (plot_name, PDFs_temp)

    ################ Draw! ################

    # Set plotting options
    x_min_str = entry_x_min.get().replace(' ', '')
    if x_min_str != '':
        x_min = float(x_min_str)
    x_max_str = entry_x_max.get().replace(' ', '')
    if x_max_str != '':
        x_max = float(x_max_str)
    
    plot1.set_xlim([x_min, x_max])
    if var_y_log.get():
        plot1.set_yscale('log')
    plot1.legend(loc='best')

    # Plot
    print('Plotting!')
    canvas.draw()

def plot_alphaN_classifier(fig, canvas, plotting_info, inputs, info):
    # Unpack information
    cuts_entries, checkboxes = inputs
    entry_x_min, entry_x_max, entry_Nbins, var_y_log = plotting_info

    ################ Setup ################

    # Clear figure
    fig.clear()

    # adding the subplot
    plot1 = fig.add_subplot(111)

    # Set labels
    plot1.set_xlabel('Log-likelyhood difference')
    plot1.set_ylabel('Normalised to global cut classifier historgram')
    plot1.set_title('Alpha-n Classifier Results')

    # Data type
    info_type = 'recon'

    colours = {
        'IBD' : ['r', 'tab:orange', 'y', 'tab:purple', 'tab:pink'],
        'alpha-n': ['b', 'tab:blue', 'tab:cyan', 'tab:gray', 'k']
    }
    linestyles = {
        'oldPDF_wo_oPs': '-',
        'newPDF_wo_oPs': '-',
        'oldPDF_wi_oPs': ':',
        'newPDF_wi_oPs': ':'
    }

    x_min_str = entry_x_min.get().replace(' ', '')
    x_max_str = entry_x_max.get().replace(' ', '')
    Nbins = entry_Nbins.get().replace(' ', '')
    if Nbins == '':
        Nbins = 20
    else:
        Nbins = int(Nbins)

    ################ plotting ################

    for process in checkboxes['process']:
        if not checkboxes['process'][process].get():
            continue
        for plot_name in checkboxes['Class']:
            if not checkboxes['Class'][plot_name].get():
                continue

            linestyle = linestyles[plot_name]
            label = process + ' ' + plot_name
            plot_name = process + '_' + plot_name

            # Upack info
            file_info, times, time_bin_lims = info[plot_name]
            file_info = file_info[info_type]

            ################ Get global cut hist, for normalising ################

            glob_cut_info = apply_global_cuts(cuts_entries['global'], file_info)
            if x_min_str != '':
                x_min = float(x_min_str)
            else:
                x_min = np.min(glob_cut_info['prompt']['Class'])
            if x_max_str != '':
                x_max = float(x_max_str)
            else:
                x_max = np.max(glob_cut_info['prompt']['Class'])
            glb_cut_hist_counts, bins = np.histogram(glob_cut_info['prompt']['Class'], bins = Nbins, range=(x_min, x_max))
            total = np.sum(glb_cut_hist_counts)
        
            ################ plotting the graph with no cuts ################

            if checkboxes['cuts']['noCUTS'].get():
                counts, _bins = np.histogram(file_info['prompt']['Class'], bins = bins)
                plot1.stairs(counts / total, bins, color=colours[process][3], linestyle = linestyle,
                             label = label + ', no cuts (N_evts = {})'.format(file_info['prompt']['E'].size))

            ################ plot with general cuts ################

            if checkboxes['cuts']['glbCUTS'].get():
                plot1.stairs(glb_cut_hist_counts / total, bins, color=colours[process][4], linestyle=linestyle,
                             label = label + ', global cuts, N_evts = {}'.format(glob_cut_info['prompt']['E'].size))

            ################ plot with PDF-specific cuts ################

            spec_cut_info = apply_PDF_specific_cuts(cuts_entries['specific'], glob_cut_info)

            for i, Class_info in enumerate(spec_cut_info):
                if checkboxes['cuts']['PDF' + str(i+1)].get():
                    counts, _bins = np.histogram(Class_info['prompt']['Class'], bins = bins)
                    plot1.stairs(counts / total, bins, color=colours[process][i], linestyle = linestyle,
                                label = label + ', PDF {}, N_evts = {}'.format(i+1, Class_info['prompt']['E'].size))

    ################ Draw! ################

    if var_y_log.get():
        plot1.set_yscale('log')
    plot1.legend(loc='best')

    # Plot
    print('Plotting!')
    canvas.draw()


####################  PRINT FUNCTIONS  ####################

def print_new_PDF(info, old_info, args):
    global PDFs

    # Setup
    f = open(args.ratdb_file)
    lines = f.readlines()
    f.close()

    # Unpack info
    old_time_bins, old_times, old_PDFs = old_info
    for process in PDFs:
        plot_name, PDFs_list = PDFs[process]
        file_info, times, time_bin_lims = info[plot_name]
        print('Printing {} to new ratdb file.'.format(plot_name))

        # Set up strings to replace in old ratdb file
        start_idx = np.argmin(np.abs(times - old_times[0]))
        end_idx = np.argmin(np.abs(times - old_times[-1])) + 1

        times = times[start_idx:end_idx]
        print('times.size =', times.size)
        PDFs_trunc = []
        for i, PDF in enumerate(PDFs_list):
            PDFs_trunc.append(PDF[start_idx:end_idx])
            print('PDFs_trunc[i].size =', PDFs_trunc[i].size)

        times_str = '['
        PDFs_str = ['[', '[', '[']
        for i in range(times.size - 1):
            times_str += str(times[i]) + ', '
            for j in range(3):
                PDFs_str[j] += str(PDFs_trunc[j][i]) + ', '
        times_str += str(times[-1]) + '],\n'
        for j in range(3):
            PDFs_str[j] += str(PDFs_trunc[j][-1]) + '],\n'

        # Replace old PDFs with new ones
        in_correct_table = False
        for i, line in enumerate(lines):
            if 'index:' in line:
                if args.index in line:
                    in_correct_table = True
                else:
                    in_correct_table = False
            elif in_correct_table:
                if 'times:' in line:
                    line = 'times: ' + times_str
                else:
                    for j in range(3):
                        string = process + str(j) + ':'
                        if string in line:
                            line = string + ' ' + PDFs_str[j]
            lines[i] = line

    # Write to brand new file
    f = open(args.new_file, 'w')
    f.writelines(lines)
    f.close()
    print('Printed to file!')


####################  GUI SETUP FUNCTIONS  ####################

def which_plot_checkboxes(Input_space):
    # Create framework
    Input_space.rowconfigure([0, 1, 2, 4], minsize=10, weight=1)
    Input_space.columnconfigure(0, minsize=10, weight=1)

    var_PDFs_or_Classifiers = tk.StringVar(Input_space, 'PDF')

    # Make title
    title_space = tk.Frame(Input_space, width=50)
    title_space.grid(row=0, column=0, sticky='nsew')
    title = tk.Label(master=title_space, text='Plot:')
    title.pack(side=tk.BOTTOM)

    # Checkbox_1 framework
    frame_1 = tk.Frame(Input_space, width=50, borderwidth=0.5, relief="ridge")
    frame_1.grid(row=1, column=0, sticky='nsew')
    frame_1.rowconfigure([0, 1], minsize=10, weight=1)
    frame_1.columnconfigure(0, minsize=10, weight=1)

    frame1_title = tk.Frame(frame_1, width=50)
    frame1_title.grid(row=0, column=0, sticky='nsew')
    tk.Radiobutton(frame1_title, text='PDFs', variable=var_PDFs_or_Classifiers, value='PDF').pack()

    frame1_options = tk.Frame(frame_1, width=50)
    frame1_options.grid(row=1, column=0, sticky='nsew')
    frame1_options.rowconfigure(0, minsize=10, weight=1)
    frame1_options.columnconfigure([0, 1, 2], minsize=10, weight=1)

    var_ratdb = tk.IntVar(value=0)
    box_ratdb = tk.Checkbutton(frame1_options, text='ratdb', variable=var_ratdb)
    box_ratdb.grid(row=0, column=0, sticky='nsew')
    var_without_oPs = tk.IntVar(value=0)
    box_without_oPs = tk.Checkbutton(frame1_options, text='without o-Ps', variable=var_without_oPs)
    box_without_oPs.grid(row=0, column=1, sticky='nsew')
    var_with_oPs = tk.IntVar(value=0)
    box_with_oPs = tk.Checkbutton(frame1_options, text='with o-Ps', variable=var_with_oPs)
    box_with_oPs.grid(row=0, column=2, sticky='nsew')

    # Checkbox_2 framework
    frame_2 = tk.Frame(Input_space, width=50, borderwidth=0.5, relief="ridge")
    frame_2.grid(row=2, column=0, sticky='nsew')
    frame_2.rowconfigure([0, 1], minsize=10, weight=1)
    frame_2.columnconfigure(0, minsize=10, weight=1)

    frame2_title = tk.Frame(frame_2, width=50)
    frame2_title.grid(row=0, column=0, sticky='nsew')
    tk.Radiobutton(frame2_title, text='Classifier Results', variable=var_PDFs_or_Classifiers, value='Class').pack()

    frame2_options = tk.Frame(frame_2, width=50)
    frame2_options.grid(row=1, column=0, sticky='nsew')
    frame2_options.rowconfigure([0, 1], minsize=10, weight=1)
    frame2_options.columnconfigure([0, 1], minsize=10, weight=1)

    var_oldPDF_wo_oPs = tk.IntVar(value=0)
    box_oldPDF_wo_oPs = tk.Checkbutton(frame2_options, text='old PDF, without o-Ps', variable=var_oldPDF_wo_oPs)
    box_oldPDF_wo_oPs.grid(row=0, column=0, sticky='nsew')
    var_oldPDF_wi_oPs = tk.IntVar(value=0)
    box_oldPDF_wi_oPs = tk.Checkbutton(frame2_options, text='old PDF, with o-Ps', variable=var_oldPDF_wi_oPs)
    box_oldPDF_wi_oPs.grid(row=0, column=1, sticky='nsew')
    var_newPDF_wo_oPs = tk.IntVar(value=0)
    box_newPDF_wo_oPs = tk.Checkbutton(frame2_options, text='new PDF, without o-Ps', variable=var_newPDF_wo_oPs)
    box_newPDF_wo_oPs.grid(row=1, column=0, sticky='nsew')
    var_newPDF_wi_oPs = tk.IntVar(value=0)
    box_newPDF_wi_oPs = tk.Checkbutton(frame2_options, text='new PDF, with o-Ps', variable=var_newPDF_wi_oPs)
    box_newPDF_wi_oPs.grid(row=1, column=1, sticky='nsew')

    #  Checkbox_3 framework
    frame_3 = tk.Frame(Input_space, width=50, borderwidth=0.5, relief="ridge")
    frame_3.grid(row=3, column=0, sticky='nsew')
    frame_3.rowconfigure(0, minsize=10, weight=1)
    frame_3.columnconfigure([0, 1], minsize=10, weight=1)

    var_alphaN = tk.IntVar(value=0)
    box_alphaN = tk.Checkbutton(frame_3, text='alpha-n', variable=var_alphaN)
    box_alphaN.grid(row=0, column=0, sticky='nsew')
    var_IDB = tk.IntVar(value=0)
    box_IBD = tk.Checkbutton(frame_3, text='IBD', variable=var_IDB)
    box_IBD.grid(row=0, column=1, sticky='nsew')

    # Checkbox_4 framework
    frame_4 = tk.Frame(Input_space, width=50, borderwidth=0.5, relief="ridge")
    frame_4.grid(row=4, column=0, sticky='nsew')
    frame_4.rowconfigure([0, 1, 2], minsize=10, weight=1)
    frame_4.columnconfigure(0, minsize=10, weight=1)

    frame4_title = tk.Frame(frame_4, width=50)
    frame4_title.grid(row=0, column=0, sticky='nsew')
    tk.Label(frame4_title, text='Cuts').pack()

    frame4_line1 = tk.Frame(frame_4, width=50)
    frame4_line1.grid(row=1, column=0, sticky='nsew')
    frame4_line1.rowconfigure(0, minsize=10, weight=1)
    frame4_line1.columnconfigure([0, 1], minsize=10, weight=1)
    var_noCUTS = tk.IntVar(value=0)
    box_noCUTS = tk.Checkbutton(frame4_line1, text='no cuts', variable=var_noCUTS)
    box_noCUTS.grid(row=0, column=0, sticky='nsew')
    var_glbCUTS = tk.IntVar(value=0)
    box_glbCUTS = tk.Checkbutton(frame4_line1, text='global cuts', variable=var_glbCUTS)
    box_glbCUTS.grid(row=0, column=1, sticky='nsew')

    frame4_line2 = tk.Frame(frame_4, width=50)
    frame4_line2.grid(row=2, column=0, sticky='nsew')
    frame4_line2.rowconfigure(0, minsize=10, weight=1)
    frame4_line2.columnconfigure([0, 1, 2], minsize=10, weight=1)
    var_PDF1 = tk.IntVar(value=0)
    box_PDF1 = tk.Checkbutton(frame4_line2, text='PDF 1', variable=var_PDF1)
    box_PDF1.grid(row=0, column=0, sticky='nsew')
    var_PDF2 = tk.IntVar(value=0)
    box_PDF2 = tk.Checkbutton(frame4_line2, text='PDF 2', variable=var_PDF2)
    box_PDF2.grid(row=0, column=1, sticky='nsew')
    var_PDF3 = tk.IntVar(value=0)
    box_PDF3 = tk.Checkbutton(frame4_line2, text='PDF 3', variable=var_PDF3)
    box_PDF3.grid(row=0, column=2, sticky='nsew')
    
    # Package outputs
    checkboxes = {}
    checkboxes['which_plot'] = var_PDFs_or_Classifiers
    checkboxes['PDF'] = {}
    checkboxes['PDF']['ratdb'] = var_ratdb; checkboxes['PDF']['without_oPs'] = var_without_oPs; checkboxes['PDF']['with_oPs'] = var_with_oPs
    checkboxes['Class'] = {}
    checkboxes['Class']['oldPDF_wo_oPs'] = var_oldPDF_wo_oPs; checkboxes['Class']['oldPDF_wi_oPs'] = var_oldPDF_wi_oPs
    checkboxes['Class']['newPDF_wo_oPs'] = var_newPDF_wo_oPs; checkboxes['Class']['newPDF_wi_oPs'] = var_newPDF_wi_oPs
    checkboxes['process'] = {}
    checkboxes['process']['alpha-n'] = var_alphaN; checkboxes['process']['IBD'] = var_IDB
    checkboxes['cuts'] = {}
    checkboxes['cuts']['noCUTS'] = var_noCUTS; checkboxes['cuts']['glbCUTS'] = var_glbCUTS
    checkboxes['cuts']['PDF1'] = var_PDF1; checkboxes['cuts']['PDF2'] = var_PDF2; checkboxes['cuts']['PDF3'] = var_PDF3

    return checkboxes

def create_input_section(Input_space):
    ################ SETUP ################
    # Create framework
    Input_space.rowconfigure([0, 1, 2, 3], minsize=50, weight=1)
    Input_space.columnconfigure(0, minsize=50, weight=1)

    # Divide space into two: global cuts, and PDF-specific cuts
    glob_title = tk.Frame(Input_space, width=50)
    glob = tk.Frame(Input_space, width=50)
    spec_title = tk.Frame(Input_space, width=50)
    spec = tk.Frame(Input_space, width=50)

    # Make titles
    label_glob_title = tk.Label(master=glob_title, text='Global Cuts:')
    label_glob_title.pack(side=tk.BOTTOM)
    label_spec_title = tk.Label(master=spec_title, text='PDF Specific Cuts (prompt E, MeV):')
    label_spec_title.pack(side=tk.BOTTOM)

    ################ GLOBAL ################
    # Create frameworks
    glob.rowconfigure([0, 1, 2, 3, 4, 5, 6], minsize=40, weight=1)
    glob.columnconfigure([0, 1, 2], minsize=30, weight=1)

    # Create widgets
    label_G_min = tk.Label(master=glob, text='Min'); label_G_max = tk.Label(master=glob, text='Max')

    label_G_promptE = tk.Label(master=glob, text='Prompt E (MeV)')
    label_G_delayedE = tk.Label(master=glob, text='Delayed E (MeV)')
    label_G_R = tk.Label(master=glob, text='Both R (m)')
    label_G_Z = tk.Label(master=glob, text='Both Z (m)')
    label_G_deltaR = tk.Label(master=glob, text='Delta R (m)')
    label_G_deltaT = tk.Label(master=glob, text='Delta T (ns)')

    entry_G_11 = tk.Entry(master=glob, width=5); entry_G_12 = tk.Entry(master=glob, width=5)
    entry_G_21 = tk.Entry(master=glob, width=5); entry_G_22 = tk.Entry(master=glob, width=5)
    entry_G_31 = tk.Entry(master=glob, width=5); entry_G_32 = tk.Entry(master=glob, width=5)
    entry_G_41 = tk.Entry(master=glob, width=5); entry_G_42 = tk.Entry(master=glob, width=5)
    entry_G_51 = tk.Entry(master=glob, width=5); entry_G_52 = tk.Entry(master=glob, width=5)
    entry_G_61 = tk.Entry(master=glob, width=5); entry_G_62 = tk.Entry(master=glob, width=5)

    # Layout
    label_G_min.grid(row=0, column=1, sticky='nsew'); label_G_max.grid(row=0, column=2, sticky='nsew')

    label_G_promptE.grid(row=1, column=0, sticky='nsew')
    label_G_delayedE.grid(row=2, column=0, sticky='nsew')
    label_G_R.grid(row=3, column=0, sticky='nsew')
    label_G_Z.grid(row=4, column=0, sticky='nsew')
    label_G_deltaR.grid(row=5, column=0, sticky='nsew')
    label_G_deltaT.grid(row=6, column=0, sticky='nsew')

    entry_G_11.grid(row=1, column=1, sticky='nsew'); entry_G_12.grid(row=1, column=2, sticky='nsew')
    entry_G_21.grid(row=2, column=1, sticky='nsew'); entry_G_22.grid(row=2, column=2, sticky='nsew')
    entry_G_31.grid(row=3, column=1, sticky='nsew'); entry_G_32.grid(row=3, column=2, sticky='nsew')
    entry_G_41.grid(row=4, column=1, sticky='nsew'); entry_G_42.grid(row=4, column=2, sticky='nsew')
    entry_G_51.grid(row=5, column=1, sticky='nsew'); entry_G_52.grid(row=5, column=2, sticky='nsew')
    entry_G_61.grid(row=6, column=1, sticky='nsew'); entry_G_62.grid(row=6, column=2, sticky='nsew')

    ################ SPECIFIC ################
    # Create frameworks
    spec.rowconfigure([0, 1, 2, 3], minsize=40, weight=1)
    spec.columnconfigure([0, 1, 2], minsize=30, weight=1)

    # Create widgets
    label_S_min = tk.Label(master=spec, text='Min'); label_S_max = tk.Label(master=spec, text='Max')

    label_S_1 = tk.Label(master=spec, text='PDF 1')
    label_S_2 = tk.Label(master=spec, text='PDF 2')
    label_S_3 = tk.Label(master=spec, text='PDF 3')

    entry_S_11 = tk.Entry(master=spec, width=5); entry_S_12 = tk.Entry(master=spec, width=5)
    entry_S_21 = tk.Entry(master=spec, width=5); entry_S_22 = tk.Entry(master=spec, width=5)
    entry_S_31 = tk.Entry(master=spec, width=5); entry_S_32 = tk.Entry(master=spec, width=5)

    # Layout
    label_S_min.grid(row=0, column=1, sticky='nsew'); label_S_max.grid(row=0, column=2, sticky='nsew')

    label_S_1.grid(row=1, column=0, sticky='nsew')
    label_S_2.grid(row=2, column=0, sticky='nsew')
    label_S_3.grid(row=3, column=0, sticky='nsew')

    entry_S_11.grid(row=1, column=1, sticky='nsew'); entry_S_12.grid(row=1, column=2, sticky='nsew')
    entry_S_21.grid(row=2, column=1, sticky='nsew'); entry_S_22.grid(row=2, column=2, sticky='nsew')
    entry_S_31.grid(row=3, column=1, sticky='nsew'); entry_S_32.grid(row=3, column=2, sticky='nsew')

    ################ OVERALL LAYOUT ################

    glob_title.grid(row=0, column=0, sticky='nsew')
    glob.grid(row=1, column=0, sticky='nsew')
    spec_title.grid(row=2, column=0, sticky='nsew')
    spec.grid(row=3, column=0, sticky='nsew')

    ################ PACKAGE ENTRIES ################

    cuts_entries = {}
    cuts_entries['global'] = {}
    cuts_entries['global']['Prompt E'] = [entry_G_11, entry_G_12]
    cuts_entries['global']['Delayed E'] = [entry_G_21, entry_G_22]
    cuts_entries['global']['Both R'] = [entry_G_31, entry_G_32]
    cuts_entries['global']['Both Z'] = [entry_G_41, entry_G_42]
    cuts_entries['global']['Delta R'] = [entry_G_51, entry_G_52]
    cuts_entries['global']['Delta T'] = [entry_G_61, entry_G_62]
    cuts_entries['specific'] = {}
    cuts_entries['specific']['PDF 1'] = [entry_S_11, entry_S_12]
    cuts_entries['specific']['PDF 2'] = [entry_S_21, entry_S_22]
    cuts_entries['specific']['PDF 3'] = [entry_S_31, entry_S_32]

    # Set default values:
    cuts_entries['global']['Prompt E'][0].insert(tk.END, '0.9')  # min (MeV)
    cuts_entries['global']['Prompt E'][1].insert(tk.END, '8.0')  # max (MeV)
    cuts_entries['global']['Delayed E'][0].insert(tk.END, '1.85')  # min (MeV)
    cuts_entries['global']['Delayed E'][1].insert(tk.END, '2.4')  # max (MeV)
    cuts_entries['global']['Both R'][1].insert(tk.END, '5.7')  # max (m)
    cuts_entries['global']['Both Z'][0].insert(tk.END, '0.85')  # min (m)
    cuts_entries['global']['Delta R'][1].insert(tk.END, '1.5')  # max (m)
    cuts_entries['global']['Delta T'][0].insert(tk.END, '400')  # min (ns)
    cuts_entries['global']['Delta T'][1].insert(tk.END, '0.8E6')  # max (ns)

    cuts_entries['specific']['PDF 1'][1].insert(tk.END, '3.5')  # max (MeV)
    cuts_entries['specific']['PDF 2'][0].insert(tk.END, '3.5')  # min (MeV)
    cuts_entries['specific']['PDF 2'][1].insert(tk.END, '5.4')  # max (MeV)
    cuts_entries['specific']['PDF 3'][0].insert(tk.END, '5.4')  # min (MeV)

    return cuts_entries
    
def make_plot_section(Input_space):
    # Creating plot frame
    plot = tk.Frame(Input_space, width=50)
    plot.rowconfigure(0, minsize=50, weight=1)
    plot.columnconfigure(0, minsize=50, weight=1)

    # creating the Tkinter canvas containing the Matplotlib figure
    fig = Figure(figsize = (5, 5), dpi = 100)
    plot1 = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master = plot)
    canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')
    canvas.draw()

    # Plotting options bar (with toolbar and other stuff)
    plotting_options_frame = tk.Frame(Input_space, width=50)
    plotting_options_frame.grid(row=1, column=0, sticky='nsew')
    plotting_options_frame.rowconfigure(0, minsize=50, weight=1)
    plotting_options_frame.columnconfigure([0, 1], minsize=50, weight=1)

    # creating the Matplotlib toolbar
    toolbar_frame = tk.Frame(plotting_options_frame, width=50)
    toolbar_frame.grid(row=0, column=1, sticky='nsew')
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)

    # Creating extra plotting options
    extra_plotting_frame = tk.Frame(plotting_options_frame, width=50)
    extra_plotting_frame.grid(row=0, column=0, sticky='nsew')
    extra_plotting_frame.rowconfigure(0, minsize=20, weight=1)
    extra_plotting_frame.columnconfigure([0, 1, 2, 3, 4, 5, 6, 7], minsize=50, weight=1)

    label_x_min = tk.Label(extra_plotting_frame, text='x min:'); label_x_min.grid(row=0, column=0, sticky='nsew')
    entry_x_min = tk.Entry(extra_plotting_frame, width=5); entry_x_min.grid(row=0, column=1, sticky='nsew')
    label_x_max = tk.Label(extra_plotting_frame, text='x max:'); label_x_max.grid(row=0, column=2, sticky='nsew')
    entry_x_max = tk.Entry(extra_plotting_frame, width=5); entry_x_max.grid(row=0, column=3, sticky='nsew')

    label_Nbins = tk.Label(extra_plotting_frame, text='N bins:'); label_Nbins.grid(row=0, column=4, sticky='nsew')
    entry_Nbins = tk.Entry(extra_plotting_frame, width=5); entry_Nbins.grid(row=0, column=5, sticky='nsew')

    label_y_log = tk.Label(extra_plotting_frame, text='y log:'); label_y_log.grid(row=0, column=6, sticky='nsew')
    var_y_log = tk.IntVar()
    box_y_log = tk.Checkbutton(extra_plotting_frame, variable=var_y_log); box_y_log.grid(row=0, column=7, sticky='nsew')

    plotting_info = (entry_x_min, entry_x_max, entry_Nbins, var_y_log)

    return plot, plotting_options_frame, fig, canvas, plotting_info


def make_GUI(info, old_info, args):
    global PDFs

    # Set up the window
    window = tk.Tk()
    window.title('Alpha_N PDFs')
    pw = ttk.PanedWindow(orient='horizontal')
    sidebar_1 = ttk.PanedWindow(pw, orient='vertical')
    sidebar_2 = ttk.PanedWindow(pw, orient='vertical')


    ###################### Plot Section ######################

    plot, plotting_options_frame, fig, canvas, plotting_info = make_plot_section(sidebar_2)


    ###################### Input Section ######################

    Input_space = tk.Frame(sidebar_1, width=50)
    Input_space.rowconfigure([0, 1], minsize=50, weight=1)
    Input_space.columnconfigure(0, minsize=50, weight=1)

    # Create Cut Input Values
    Cut_space = tk.Frame(Input_space, width=50)
    cuts_entries = create_input_section(Cut_space)

    # Create plotting checkboxes
    Checkbox_space = tk.Frame(Input_space, width=50, borderwidth=0.5, relief="raised")
    checkboxes = which_plot_checkboxes(Checkbox_space)

    # Arrange sections, and package output together
    Cut_space.grid(row=0, column=0, sticky='nsew')
    Checkbox_space.grid(row=1, column=0, sticky='nsew')

    ###################### Buttons Section ######################
    Buttons_space = tk.Frame(sidebar_1, width=20)
    Buttons_space.rowconfigure(0, minsize=10, weight=1)
    Buttons_space.columnconfigure([0, 1], minsize=10, weight=1)

    # Create Plot button
    inputs = (cuts_entries, checkboxes)
    plot_button = tk.Button(Buttons_space, text='Plot',
                            command=lambda: plotting(fig, canvas, plotting_info, inputs, info, old_info))
    plot_button.grid(row=0, column=0, sticky='nsew')

    # Create Print button
    print_button = tk.Button(Buttons_space, text='Print PDFs',
                            command=lambda: print_new_PDF(info, old_info, args))
    print_button.grid(row=0, column=1, sticky='nsew')

    ###################### Finalising Section ######################

    # add the two sidebars to the main window
    pw.add(sidebar_1)
    pw.add(sidebar_2)

    # add the two halves of the left sidebar to the sidebar
    sidebar_1.add(Input_space)
    sidebar_1.add(Buttons_space)

    # add the two halves of the right sidebar to the sidebar
    sidebar_2.add(plot)
    sidebar_2.add(plotting_options_frame)

    # add the paned window to the root
    pw.pack(fill='both', expand=True)

    window.mainloop()


def main():
    args = argparser()

    info_file_names = ['IBD_oldPDF_wo_oPs', 'alpha-n_oldPDF_wo_oPs', 'IBD_oldPDF_wi_oPs', 'alpha-n_oldPDF_wi_oPs',
                       'IBD_newPDF_wo_oPs', 'alpha-n_newPDF_wo_oPs', 'IBD_newPDF_wi_oPs', 'alpha-n_newPDF_wi_oPs']
    info = {}
    for i, info_file in enumerate(args.info_files):
        if os.path.exists(info_file):
            # Read in info from text file to appropriate data structure
            file_info, times, time_bins = read_info(info_file)

            # Transform into prompt-delayed event structure
            re_info = create_prompt_delay(file_info)

            # Package
            info[info_file_names[i]] = (re_info, times, time_bins)

    # Read in ratdb/old_PDF file
    old_time_bins, old_times, old_PDFs = read_ratdb_file(args.ratdb_file, args.index)
    old_info = (old_time_bins, old_times, old_PDFs)

    # Make GUI and add plot
    make_GUI(info, old_info, args)


if __name__ == '__main__':
    main()