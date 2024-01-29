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
                        # default='macros/labppo_2p2_scintillator/flat/alphaN_13C_flat.mac')
                        # default='macros/labppo_2p2_scintillator/flat/IBD_flat.mac')
                        default='macros/labppo_2p2_scintillator/AmBe.mac')

    parser.add_argument('--sim_repo', '-sr', type=str, dest='sim_repo',
                        # default='/mnt/lustre/scratch/epp/jp643/antinu/AmBe/RAT7.0.15/sims/', help='Folder to save intial root files from simulations in.')
                        default='/mnt/lustre/scratch/epp/jp643/antinu/Analysis_data/AmBe/RAT7.0.15/', help='Folder to save intial root files from simulations in.')
                        # default='/mnt/lustre/scratch/epp/jp643/antinu/MC_data/AmBe/RAT7.0.15/', help='Folder to save intial root files from simulations in.')
    parser.add_argument('--info_repo', '-ir', type=str, dest='info_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/AmBe/RAT7.0.15/info/', help='Folder to save info text files in.')
    parser.add_argument('--hist_repo', '-hr', type=str, dest='hist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/AmBe/RAT7.0.15_Nhit/hists/', help='Folder to save hist root files in.')
    parser.add_argument('--tothist_repo', '-tr', type=str, dest='tothist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/AmBe/RAT7.0.15_Nhit/tothists/', help='Folder to save combined total hist root files in.')
    
    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=200000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=500, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--max_jobs', '-mx', type=int, dest='max_jobs',
                        default=30, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---step', '-s', type=str, dest='step', required=True, choices=['sim', 'resim', 'info', 'hist', 'combi'],
                        help='which step of the process is it in?')
    parser.add_argument('---start_fileNum', '-sn', type=int, dest='start_fileNum', default=0,
                        help='Which number the files are numbered from (for splitting simulations into multiple files for example)')
    
    parser.add_argument('---runNum_list', '-rl', type=str, dest='runNum_list', help='Text file with run numbers to simulate (one per line)',
                        # default='')
                        default='run_lists/AmBe.txt')
    parser.add_argument('---use_all_files_in_dir', '-A', type=bool, dest='use_all_files_in_dir',
                        default=True, help='For [info] or [hist] modes. Instead of using the number of events to work out which files to get info from in sim directory, just use all the files in there.')
    parser.add_argument('---flat', '-f', type=bool, dest='flat',
                        default=False, help='True if you want to also produce the same histograms, but where the prompt E spectra have been flattened.')
    parser.add_argument('---is_data', '-iD', type=bool, dest='is_data',
                        default=True, help='For energy correction: True for data, False for MC.')
    
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                        default=False, help='print and save extra info')

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
            new_line = '/rat/run/start\n'  # Don't say how many events. Will be set in runtime argument
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
    new_macro_address = jobScript_repo + 'macro_' + filename_format(args.macro) + '.mac'
    commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(args.macro) + '.txt'
    new_job_address = jobScript_repo + job_str_map('sims_', args.macro) + '.job'
    output_logFile_address = logFile_repo + 'log_sims_' + filename_format(args.macro) + '.txt'

    return save_sims_folder, logFile_repo, new_macro_address, commandList_address, new_job_address, output_logFile_address

def make_command_list(save_sims_folder, logFile_repo, new_macro_address, n_evts, run_numbers, args):
    '''Create list of commands to be printed to file, and used by job script.'''

    commands = []
    for run in run_numbers:
        for i in range(len(n_evts)):
            # Create all the commands to run the macro
            outroot_address = save_sims_folder + 'simOut_' + filename_format(args.macro)
            log_file_address = logFile_repo + 'ratLog_' + filename_format(args.macro)

            if run:
                outroot_address += '_' + str(run)
                log_file_address += '_' + str(run)

            outroot_address += '_' + str(args.start_fileNum + i) + '.root'
            log_file_address += '_' + str(args.start_fileNum + i) + '.log'

            macro_command = 'rat ' + new_macro_address + ' -N ' + str(n_evts[i]) + ' -o ' + outroot_address + ' -l ' + log_file_address

            if run:
                macro_command += ' -n ' + str(run)  # Choose the run number (can also choose sub-run number)
            else:
                macro_command += ' -P'  # Airplane mode (no remote db access, since unnecessary without a specific run)
                
            if args.verbose:
                macro_command += ' -vv'
            
            commands.append(macro_command + '\n')
    return commands

def runSims(args):
    '''Runs simulations based in input information'''

    print('Running runSims().')
    save_sims_folder, logFile_repo, new_macro_address, commandList_address, new_job_address, output_logFile_address = make_addresses(args)

    # Read in example macro and job script + info
    repo_address = getRepoAddress()
    with open(repo_address + args.macro, 'r') as f:
        example_macro = f.readlines()

    example_jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(example_jobScript_address, 'r') as f:
        example_jobScript = f.readlines()

    # Read in run number info, if provided
    if args.runNum_list != '':
        run_numbers = []
        with open(repo_address + args.runNum_list, 'r') as f:
            for num in f.readlines():
                run_numbers.append(int(num))
    else:
        run_numbers = [False]
    print()

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
    
    ### MAKE MACROS AND JOB SCRIPTS TO RUN THE SIMULATIONS ###
    print('Creating macros and job scripts...')
    
    # Create macro
    makeMacro(new_macro_address, example_macro)

    # Create command list file to call macro repeatedly
    commands = make_command_list(save_sims_folder, logFile_repo, new_macro_address, n_evts, run_numbers, args)
    with open(commandList_address, 'w') as f:
        command_list = ''.join(commands)
        f.write(command_list)

    # Create the job script to run all these commands in a file
    job_address = makeJobArrayScript(new_job_address, output_logFile_address, example_jobScript, commandList_address)

    ### RUN JOB SCRIPTS ###
    print('Submitting job array...')
    command = 'qsub -l m_mem_free=4G -t 1-' + str(len(n_evts) * len(run_numbers)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True 

### RERUN FUNCTIONS ###

def check_failed(log_address):
    '''Check end of log file to see if simulation completed successfully. Return True if it did not.'''

    # If log file doesn't exist, it definitely failed.
    if not os.path.exists(log_address):
        return True

    with open(log_address, 'r') as f:
        lines = f.readlines()

    # Lines that should always be present towards the end of a successfully completed rat logfile.
    checks = {
        'ProcBlock::~ProcBlock Processor usage statistics': False,
        'ConditionalProcBlock::~ConditionalProcBlock Processor usage statistics': False,
        'Gsim::~Gsim Event simulation statistics': False,
        'GLG4PrimaryGeneratorAction::~GLG4PrimaryGeneratorAction: Deleting generators.': False,
        'Graphics systems deleted.': False,
        'Visualization Manager deleting...': False
    }

    for line in lines:
        stripped_line = line.strip()
        if stripped_line in checks:
            checks[stripped_line] = True

    lines_present = True
    for key in checks:
        lines_present = lines_present and checks[key]
    
    return not lines_present

def reSim(args):
    '''Check which simulations failed, and rerun them'''

    print('Running reSim().')
    _, _, _, commandList_address, new_job_address, _ = make_addresses(args)

    # checkJobsDone('sims_', args.macro, 10, args.verbose)

    # Read in command list file, to get log file addresses
    with open(commandList_address, 'r') as f:
        commands = f.readlines()

    # Get log file, and output sim file addresses
    rerun_info = []
    rerun_numbers = ''
    for i, command in enumerate(commands):
        if ' -o ' not in command or ' -l ' not in command:
            print('No outRoot file or log file defined in line: "{}"',format(command))
        else:
            arguments = command.split(' ')
            outRoot_file = ''
            log_file = ''
            for j, argument in enumerate(arguments):
                if argument == '-o':
                    outRoot_file = arguments[j + 1]
                elif argument == '-l':
                    log_file = arguments[j + 1]

            if outRoot_file != '' and log_file != '':
                # Check if simulation failed
                failed = check_failed(log_file)

                if failed:
                    rerun_info.append((i, outRoot_file, log_file))
                    rerun_numbers += ', ' + str(i)
    if rerun_numbers[:2] == ', ':
        rerun_numbers = rerun_numbers[2:]
    
    if len(rerun_info) > 0:
        print('{} simulation(s) failed: {}'.format(len(rerun_info), rerun_numbers))
        confirmation = input('Delete outputs, and rerun them? Answer with Y/!Y. ANSWER: ')
        if confirmation in ('Y', 'y', True):
            start_idx = -1
            end_idx = -1
            for i, outRoot_file, log_file in rerun_info:
                # Delete outRoot and log files
                if os.path.exists(outRoot_file):
                    if args.verbose:
                        print('Deleting: {}'.format(outRoot_file))
                    command = 'rm ' + outRoot_file
                    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished
                if os.path.exists(log_file):
                    if args.verbose:
                        print('Deleting: {}'.format(log_file))
                    command = 'rm ' + log_file
                    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

                # Try to group subjobs into the same command, if their numbers are ajacent
                if start_idx == -1:
                    start_idx = i
                    end_idx = i
                if i == end_idx + 1:
                    end_idx += 1
                else:
                    # Rerun simulation
                    command = 'qsub -l m_mem_free=4G -t ' + str(start_idx+1) + '-' + str(end_idx+1) + ' -tc ' + str(args.max_jobs) + ' ' + new_job_address
                    if args.verbose:
                        print('Running command: ', command)
                    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

                    start_idx = i
                    end_idx = i
            
            if start_idx != end_idx:
                # Rerun simulation (last group of subjobs)
                command = 'qsub -l m_mem_free=4G -t ' + str(start_idx+1) + '-' + str(end_idx+1) + ' -tc ' + str(args.max_jobs) + ' ' + new_job_address
                if args.verbose:
                    print('Running command: ', command)
                subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished
    else:
        print('No simulations failed.')

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

    command = command_base + output_file + ' ' + str(int(args.start_fileNum) * int(args.nevts_persim)) + ' ' + str(int(args.verbose))

    if args.use_all_files_in_dir:
        for filename in os.listdir(args.sim_repo):
            file_address = os.path.join(args.sim_repo, filename)
            if os.path.isfile(file_address):
                if file_address[-5:] == '.root':
                    command += ' ' + file_address
    else:
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


def getHists(args):
    '''Get info from simulaltion output and print everything to text file'''
    print('Running getHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    ### Create file and folder addresses needed by other functions ###

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_hist_folder = checkRepo(args.hist_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_hist_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # Log file folder
    logFile_repo = save_hist_folder + 'log/'
    logFile_repo = checkRepo(logFile_repo, args.verbose)

    # New addresses
    new_job_address = jobScript_repo + 'hist_job_' + filename_format(args.macro) + '.job'
    output_logFile_address = logFile_repo + 'log_hist_' + filename_format(args.macro) + '.txt'
    commandList_address = jobScript_repo + 'hist_commandList.txt'

    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    command_base = repo_address + 'scripts/make_hists.exe '

    file_addresses = []
    if args.use_all_files_in_dir:
        out_address_start = save_hist_folder + 'hist'
        for filename in os.listdir(args.sim_repo):
            file_address = os.path.join(args.sim_repo, filename)
            if os.path.isfile(file_address):
                if file_address[-5:] == '.root':
                    file_addresses.append(' ' + file_address)
    else:
        # How simulations were split up
        n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
        out_address_start = save_hist_folder + 'hist_' + filename_format(args.macro)
        for i in range(len(n_evts)):
            file_addresses.append(' ' + save_sims_folder + 'simOut_' + filename_format(args.macro) + '_' + str(args.start_fileNum + i) + '.root')
    
    if args.flat:
        file_addresses_str = ''
        for file_add in file_addresses:
            file_addresses_str += file_add
        command = command_base + out_address_start + '.root ' + out_address_start + '.txt ' + str(int(args.flat)) + ' ' + str(int(args.verbose)) + file_addresses_str + '\n'

        # Get example job script
        example_jobScript_address = repo_address + 'job_scripts/jobSingle.job'
        with open(example_jobScript_address, "r") as f:
            example_jobScript = f.readlines()
        makeJobSingleScript(new_job_address, output_logFile_address, example_jobScript, command)

        # Wait until previous jobs are done
        checkJobsDone('sims_', args.macro, 10, args.verbose)

        # Make job submission command
        command = 'qsub -l m_mem_free=4G ' + new_job_address

    else:
        commands = []
        for i in range(len(file_addresses)):
            outRoot_address = out_address_start + '_' + str(i) + '.root '
            outText_address = out_address_start + '_' + str(i) + '.txt '
            commands.append(command_base + outRoot_address + outText_address + str(int(args.flat)) + ' ' + str(int(args.is_data)) + ' ' + str(int(args.verbose)) + file_addresses[i] + '\n')
        
        # Create the job script to run all these commands in a file
        with open(commandList_address, 'w') as f:
            command_list = ''.join(commands)
            f.write(command_list)

        # Create the job script to run all these commands in a file
        example_jobScript_address = repo_address + 'job_scripts/jobArray.job'
        with open(example_jobScript_address, 'r') as f:
            example_jobScript = f.readlines()
        job_address = makeJobArrayScript(new_job_address, output_logFile_address, example_jobScript, commandList_address)

        command = 'qsub -t 1-' + str(len(commands)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 


    ### RUN JOB SCRIPTS ###
    print('Submitting job array...')
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True

def combiHists(args):
    '''Combine split histograms together'''
    print('Running combiHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_splithists_folder = checkRepo(args.hist_repo, args.verbose)
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)

    # Get example job script
    example_jobScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(example_jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # Making job scripts
    command = 'hadd ' + save_tothists_folder + 'tot_hists.root ' + save_splithists_folder + 'hist_*.root'

    # Make job script
    new_job_address = save_tothists_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, args.verbose)
    new_job_address += job_str_map('info_', args.macro) + '.job'

    output_logFile_address = save_tothists_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, args.verbose)
    output_logFile_address +=  'log_info_' + filename_format(args.macro) + '.txt'
    makeJobSingleScript(new_job_address, output_logFile_address, example_jobScript, command)
    
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
        'resim': reSim,
        'info': getInfo,
        'hist': getHists,
        'combi': combiHists
    }

    result = work_modes[args.step](args)

if __name__ == '__main__':
    main()
