import numpy as np
import argparse
import subprocess
import os
import time

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--particle', '-p', type=str, dest='particle',
                        default='o-Ps', choices=['o-Ps', 'e+', 'e-', 'gamma', 'alpha'], help='Which particle to simulate')
    parser.add_argument('--energies', '-e', type=str, dest='energies',
                        default=[1.0], help='List of particle energies to simulate (MeV), Input as [1.0,2.0,3.0] for ex')

    parser.add_argument('--macro_repo', '-mr', type=str, dest='macro_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/macros/', help='Folder to save Region-selected root files in.')
    parser.add_argument('--sim_repo', '-sr', type=str, dest='sim_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/sims/', help='Folder to save intial root files from simulations in.')
    parser.add_argument('--pdf_repo', '-pr', type=str, dest='pdf_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/PDFs/', help='Folder to save PDF text files in.')

    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=10000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=1000, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---step', '-s', type=str, dest='step', required=True, choices=['sim', 'pdf', 'sim-pdf'],
                        help='which step of the process is it in?')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                    default=True, help='print and save extra info')

    args = parser.parse_args()
    return args


### MISCELLANEOUS FUNCTIONS ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('scripts/runAMMELIE.py')]
    if repo_address == '':
        firt_char = None
    else:
        firt_char = repo_address[0]
    if firt_char != '/':
        repo_address = os.getcwd() + '/' + repo_address
    return repo_address

def checkRepo(repo_address, verbose=False):
    '''Check format of repo address is right to be used here. Also if repo does not exist, create it.'''

    # Format address
    new_address = repo_address
    if new_address == '.':
        new_address = ''
    elif new_address != '':
        if (new_address[int(len(new_address) - 1)] != '/'):
            new_address += '/'
        
        # If directory does not exists, make it
        if not os.path.isdir(new_address):
            os.mkdir(new_address)
            if verbose:
                print('Created new directory: ', new_address)

    return new_address

def getNevtsPerMacro(nevts_total, nevts_persim):
    '''Create array with number of events to simulate per macro'''

    n_macros = nevts_total // nevts_persim
    remainder = nevts_total % nevts_persim
    if remainder == 0:
        n_evts = np.an_array = np.full(n_macros, nevts_persim)
    else:
        n_evts = np.an_array = np.full(n_macros + 1, nevts_persim)
        n_evts[n_macros] = remainder

    return n_evts

def filename_format(particle, energy):
    '''returns string with simulation info to use in filenames'''
    return 'Positronium_' + particle + '_' + str(energy)

def job_str_map(jobName_str, particle, energy):
    '''Create code string to put at the start of the job file name, so that it can
    be recognised in the job list (the job list only displays the first 10 characters).'''

    map = {
        'job_type': {
            'sims_': 'S',
            'pdf_': 'P'
        },
        'particle': {
            'o-Ps': 'oP',
            'e+': 'e+',
            'e-': 'e-',
            'gamma': 'gm',
            'alpha': 'al'
        }
    }

    return map['job_type'][jobName_str] + '_' + map['particle'][particle] + '_' + str(energy)

def makeJobArrayScript(jobName_str, example_jobScript, overall_folder, commandList_address, particle, energy, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, particle, energy) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(particle, energy) + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif '/Address/CommandList.txt' in line:
            new_line = line.replace('/Address/CommandList.txt', commandList_address, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def makeJobSingleScript(jobName_str, example_jobScript, overall_folder, commands, particle, energy, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, particle, energy) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(particle, energy) + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif 'your commands' in line:
            new_line = line.replace('your commands', commands, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def checkJobsDone(jobName_str, particle, energies, wait_time, verbose):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    # Turns out the name of the job is only the 10 first characters of the job file name

    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE, shell=True).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
                for energy in energies:
                    map_str = job_str_map(jobName_str, particle, energy)
                    if len(map_str) > 10:
                        map_str = map_str[:9]
                    if map_str in line:
                        running = True
                        if verbose:
                            print('Waiting for jobs to finish...')
                        break
        time.sleep(wait_time)

    return True


### SIMS FUNCTIONS ###

def makeMacros(particle, energy, example_macro, save_macro_folder, save_sims_folder, n_evts, idx):
    '''Make and save macro according to provided parameters'''

    # AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root
    new_macro_address = save_macro_folder + 'macro_' + filename_format(particle, energy) + '_' + str(idx) + '.mac'
    output_address = save_sims_folder + 'simOut_' + filename_format(particle, energy) + '_' + str(idx) + '.root'

    new_macro = []
    for line in example_macro:
        # Replace placeholders in macro
        if ('/rat/db/set POSITRONIUM formfrac 1.0' in line) and (particle != 'o-Ps'):
            new_line = line.replace('1.0', '0.0')
        elif ('/rat/tracking/omit e-' in line) and (particle == 'e-'):
            new_line = line.replace('/rat/tracking/omit e-', '')
        elif '/rat/procset file "oPs_output.root"' in line:
            new_line = line.replace('oPs_output.root', output_address)
        elif '/generator/vtx/set e+ 0 0 0 1.0' in line:
            new_line = line.replace('e+', particle)
            new_line = line.replace('1.0', str(energy))
        elif '/rat/run/start 1000' in line:
            new_line = line.replace('1000', str(n_evts))
        else:
            new_line = line

        new_macro.append(new_line)
    
    # Create new macro file
    with open(new_macro_address, "w") as f:
        new_macro = "".join(new_macro)
        f.write(new_macro)

    return new_macro_address

def runSims(args):
    '''Runs simulations based in input information'''
    print('Running runSims().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    example_macro_address = repo_address + 'macros/Sim_oPs.mac'
    with open(example_macro_address, "r") as f:
        example_macro = f.readlines()

    example_jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(example_jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_macro_folder = checkRepo(args.macro_repo, args.verbose)
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_macro_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
    
    ### MAKE MACROS AND JOB SCRIPTS TO RUN THE SIMULATIONS ###
    print('Creating macros and job scripts...')
    job_addresses = []
    energies = map(float, args.energies.strip('[]').split(','))
    for energy in energies:
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(args.particle, energy) + '.txt'
        commandList_file = open(commandList_address, 'w')
        for i in range(len(n_evts)):
            # Create all the macros
            macro_address = makeMacros(args.particle, energy, example_macro, save_macro_folder, save_sims_folder, n_evts[i], i)
            log_file_address = save_macro_folder + 'log_files/ratLog_' + filename_format(args.particle, energy) + '.log'
            macro_command = 'rat ' + macro_address + ' -l ' + log_file_address
            commandList_file.write(macro_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('sims_', example_jobScript, save_macro_folder, commandList_address, args.particle, energy, args.verbose)
        job_addresses.append(new_job_address)


    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub -t 1-' + str(len(n_evts)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### Analysis functions ###

def getPDF(args):
    '''Compute stats'''
    print('Running getPDF().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_pdf_folder = checkRepo(args.pdf_repo, args.verbose)

    # Get example job script
    example_jobScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(example_jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # How simulations were split up
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)

    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    jobScript_addresses = []
    command_base = repo_address + 'scripts/getPDF.exe '
    energies = map(float, args.energies.strip('[]').split(','))
    for energy in energies:
        output_file = save_pdf_folder + 'PDFs_' + filename_format(args.particle, energy) + '.txt'

        command = command_base + output_file
        if args.particle == 'o-Ps':
            command += ' ' + str(1)
        else:
            command += ' ' + str(0)
        command += ' ' + str(int(args.verbose))

        sim_file_format = save_sims_folder + 'simOut_' + filename_format(args.particle, energy)
        for i in range(len(n_evts)):
            command += ' ' + sim_file_format + '_' + str(i) + '.root'

        new_job_address = makeJobSingleScript('pdf_', example_jobScript, save_sims_folder, command, args.particle, energy, args.verbose)
        jobScript_addresses.append(new_job_address)

    # Wait until previous jobs are done
    checkJobsDone('sims_', args.particle, energies, 10, args.verbose)

    ### RUN JOB SCRIPTS ###
    print('Submitting job(s)...')
    for job_address in jobScript_addresses:
        # For higher energies, the jobs need more memory, otherwise they get killed
        command = 'qsub -l m_mem_free=4G ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### COMBINATIONS ###

# def runSim_getPDF(args):
#     result = runSims(args)
#     return result and getPDF(args)

### MAIN ###

def main():
    # read in argument
    args = argparser()

    work_modes = {
        'sim': runSims,
        'pdf': getPDF,
        #'sim-pdf': runSim_getPDF
    }

    result = work_modes[args.step](args)

if __name__ == '__main__':
    main()