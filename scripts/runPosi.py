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

    parser.add_argument('--macro', '-m', type=str, dest='macro', help='Which macro to base simulation macros off of.',
                        default='macros/2p2labppo/ReactorIBD.mac')
                        # default='macros/2p2labppo/alphaN_13C.mac')

    parser.add_argument('--sim_repo', '-sr', type=str, dest='sim_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/labppo_2p2_scintillator/makePDFs/sims/', help='Folder to save intial root files from simulations in.')
    parser.add_argument('--info_repo', '-ir', type=str, dest='info_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Positronium/labppo_2p2_scintillator/makePDFs/info/', help='Folder to save info text files in.')
    
    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=10000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=1000, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---step', '-s', type=str, dest='step', required=True, choices=['sim', 'info'],
                        help='which step of the process is it in?')
    parser.add_argument('---start_fileNum', '-sn', type=int, dest='start_fileNum', default=0,
                        help='Which number the files are numbered from (for splitting simulations into multiple files for example)')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                    default=True, help='print and save extra info')

    args = parser.parse_args()
    return args


### MISCELLANEOUS FUNCTIONS ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('scripts/runPosi.py')]
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

def filename_format(macro_name):
    '''returns string with simulation info to use in filenames'''
    return macro_name.split('/')[-1].replace('.mac', '')

def job_str_map(jobName_str, macro_name):
    '''Create code string to put at the start of the job file name, so that it can
    be recognised in the job list (the job list only displays the first 10 characters).'''

    map = {
        'sims_': 'S',
        'info_': 'I'
    }

    return map[jobName_str] + '_' + filename_format(macro_name)[:9]

def makeJobArrayScript(new_job_address, output_logFile_address, example_jobScript, commandList_address):
    '''Create job script to run array of rat macros'''

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

def makeJobSingleScript(new_job_address, output_logFile_address, example_jobScript, commands):
    '''Create job script to run array of rat macros'''

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

def checkJobsDone(jobName_str, macro_name, wait_time, verbose):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    # Turns out the name of the job is only the 10 first characters of the job file name

    map_str = job_str_map(jobName_str, macro_name)
    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE, shell=True).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
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

def makeMacro(new_macro_address, example_macro):
    '''Make and save macro according to provided parameters'''

    new_macro = []
    for line in example_macro:
        if '/rat/procset file' in line:
            new_line = '#' + line  # Will set with run-time argument
        elif '/rat/run/start' in line:
            new_line = '#' + line  # Will set with run-time argument
        else:
            new_line = line

        new_macro.append(new_line)
    
    # Create new macro file
    with open(new_macro_address, 'w') as f:
        new_macro = ''.join(new_macro)
        f.write(new_macro)

def make_addresses(args):
    '''Create file and folder addresses needed by other functions.'''

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_sims_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # Log file folder
    logFile_repo = save_sims_folder + 'log/'
    logFile_repo = checkRepo(logFile_repo, args.verbose)

    # New addresses
    new_macro_address = save_sims_folder + 'macro_' + filename_format(args.macro) + '.mac'
    commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(args.macro) + '.txt'
    new_job_address = jobScript_repo + job_str_map('sims_', args.macro) + '.job'
    output_logFile_address = logFile_repo + 'log_sims_' + filename_format(args.macro) + '.txt'

    return save_sims_folder, new_macro_address, commandList_address, new_job_address, output_logFile_address

def make_command_list(save_sims_folder, new_macro_address, n_evts, args):
    '''Create list of commands to be printed to file, and used by job script.'''

    commands = []
    for i in range(len(n_evts)):
        # Create all the commands to run the macro
        outroot_address = save_sims_folder + 'simOut_' + filename_format(args.macro) + '_' + str(args.start_fileNum + i) + '.root'
        log_file_address = save_sims_folder + 'log_files/ratLog_' + filename_format(args.macro) + '_' + str(args.start_fileNum + i) + '.log'
        macro_command = 'rat -P ' + new_macro_address + ' -N ' + str(n_evts[i]) + ' -o ' + outroot_address + ' -l ' + log_file_address
        if args.verbose:
            macro_command += ' -vv'
        commands.append(macro_command)
    return commands

def runSims(args):
    '''Runs simulations based in input information'''

    print('Running runSims().')
    save_sims_folder, new_macro_address, commandList_address, new_job_address, output_logFile_address = make_addresses(args)

    # Read in example macro and job script + info
    repo_address = getRepoAddress()
    with open(repo_address + args.macro, 'r') as f:
        example_macro = f.readlines()

    example_jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(example_jobScript_address, 'r') as f:
        example_jobScript = f.readlines()

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
    
    ### MAKE MACROS AND JOB SCRIPTS TO RUN THE SIMULATIONS ###
    print('Creating macros and job scripts...')
    
    # Create macro
    makeMacro(new_macro_address, example_macro)

    # Create command list file to call macro repeatedly
    commands = make_command_list(save_sims_folder, new_macro_address, n_evts, args)
    with open(commandList_address, 'w') as f:
        command_list = ''.join(commands)
        f.write(command_list)

    # Create the job script to run all these commands in a file
    job_address = makeJobArrayScript(new_job_address, output_logFile_address, example_jobScript, commandList_address)

    ### RUN JOB SCRIPTS ###
    print('Submitting job array...')
    command = 'qsub -t 1-' + str(len(n_evts)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True 

### Analysis functions ###

def getInfo(args):
    '''Get info from simulaltion output and print everything to text file'''
    print('Running getInfo().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_info_folder = checkRepo(args.info_repo, args.verbose)

    # Get example job script
    example_jobScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(example_jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # How simulations were split up
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)

    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    command_base = repo_address + 'scripts/get_evt_info.exe '
    output_file = save_info_folder + 'info_' + filename_format(args.macro) + '_' + str(args.start_fileNum) + '.txt'

    command = command_base + output_file + ' ' + str(args.start_fileNum) + ' ' + str(int(args.verbose))

    sim_file_format = save_sims_folder + 'simOut_' + filename_format(args.macro)
    for i in range(len(n_evts)):
        command += ' ' + sim_file_format + '_' + str(args.start_fileNum + i) + '.root'

    # Make job script
    new_job_address = save_info_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, args.verbose)
    new_job_address += job_str_map('info_', args.macro) + '.job'

    output_logFile_address = save_info_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, args.verbose)
    output_logFile_address +=  'log_info_' + filename_format(args.macro) + '.txt'
    makeJobSingleScript(new_job_address, output_logFile_address, example_jobScript, command)

    # Wait until previous jobs are done
    checkJobsDone('sims_', args.macro, 10, args.verbose)

    ### RUN JOB SCRIPTS ###
    print('Submitting job...')
    # For higher energies, the jobs need more memory, otherwise they get killed
    command = 'qsub -l m_mem_free=4G ' + new_job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True

### MAIN ###

def main():
    # read in argument
    args = argparser()

    work_modes = {
        'sim': runSims,
        'info': getInfo
    }

    result = work_modes[args.step](args)

if __name__ == '__main__':
    main()
