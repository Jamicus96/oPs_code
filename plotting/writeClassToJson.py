import json
import argparse
import os
import numpy as np


def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--class_repo', '-pr', type=str, dest='class_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/Classifications/', help='Folder to save Classifier result text files in.')
    parser.add_argument('--json_repo', '-jsr', type=str, dest='json_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/Classifications/',
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
    for file_address, particle, particle_energy in info:
        file = open(file_address, 'r')
        lines = file.readlines()
        file.close()

        if particle not in table:
            table[particle] = {}
        table[particle][particle_energy]= {}
        table[particle][particle_energy]['Classifier_results'] = []
        table[particle][particle_energy]['delays'] = []
        table[particle][particle_energy]['nhits'] = []
        table[particle][particle_energy]['recon_event_energies'] = []
        for line in lines:
            stats = line.split(' ')
            table[particle][particle_energy]['Classifier_results'].append(float(stats[0]))
            table[particle][particle_energy]['delays'].append(float(stats[1]))
            table[particle][particle_energy]['nhits'].append(float(stats[2]))
            table[particle][particle_energy]['recon_event_energies'].append(float(stats[3]))

    with open(json_repo + 'Classifier_stats.json', 'w') as f:
        json.dump(table, f)


### MAIN ###

def main():
    # read in argument
    args = argparser()

    # List files in repo
    files = listDir(args.class_repo)

    # Get all relevent info from file name
    info = getNameInfo(args.class_repo, files)

    # Write PDFs to json table
    writeToJson(args.json_repo, info)

if __name__ == '__main__':
    main()