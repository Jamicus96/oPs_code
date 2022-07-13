import json
import argparse
import os
import numpy as np


def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--pdf_repo', '-pr', type=str, dest='pdf_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/PDFs/', help='Folder to save PDF text files in.')
    parser.add_argument('--json_repo', '-jsr', type=str, dest='json_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/',
                        help='Folder to json file with final stats in.')

    args = parser.parse_args()
    return args

def listDir(dir):
    '''List all the files in a directory'''
    init_list = os.listdir(dir)
    file_list = []
    for item in init_list:
        if os.path.isfile(dir + item):
            file_list.append(item)
    return file_list

def getNameInfo(dir, files):

    info = []
    for file in files:
        file_info = file.replace('.txt', '').split('_')
        particle = file_info[2]
        energy = file_info[3]
        if particle == '' or energy == '' or 'Hists_' in file:
            continue 
        info.append((dir + file, particle, energy))
    return info

def writeToJson(json_repo, info):

    table = {}
    for file_address, particle, energy in info:
        file = open(file_address, 'r')
        lines = file.readlines()
        file.close()

        times = list(np.array(lines[0].replace('[', '', 1).replace(']', '', 1).split(', ')).astype(float))
        MC_probs = list(np.array(lines[1].replace('[', '', 1).replace(']', '', 1).split(', ')).astype(float))
        Fitted_probs = list(np.array(lines[2].replace('[', '', 1).replace(']', '', 1).split(', ')).astype(float))

        if particle not in table:
            table[particle] = {}
        if 'times' not in table:  # Assuming times are the same for all files/PDFs
            table['times'] = times

        table[particle][energy] = {}
        table[particle][energy]['PDF_MC'] = MC_probs
        table[particle][energy]['PDF_Fitted'] = Fitted_probs
    with open(json_repo + 'PDF_stats.json', 'w') as f:
        json.dump(table, f)


### MAIN ###

def main():
    # read in argument
    args = argparser()

    # List files in repo
    files = listDir(args.pdf_repo)

    # Get all relevent info from file name
    info = getNameInfo(args.pdf_repo, files)

    # Write PDFs to json table
    writeToJson(args.json_repo, info)

if __name__ == '__main__':
    main()