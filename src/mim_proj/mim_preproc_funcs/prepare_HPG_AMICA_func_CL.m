function [txt] = prepare_HPG_AMICA_func_CL(EEG,fileNameNoExt,amicaOutputFolder_local,amicaOutputFolder_unix,avgRefPCAReduction,emailStr)
%prepare_HPG_AMICA_func help
%
%EEG is the EEG struct you want to save on M drive and run AMICA on via the
%hipergator
%
%fileNameNoExt is a string for what you want to save the file as (without
%.set extension)
%
%amicaOutputFolder is a string for where you want to store the eeg data set
%and accompanying parameter files. You just need to supply the general 
%location like 'M:\share\MindInMotion\SUBJ####\EEG\AMICA' or
%fullfile(MIMDataFolder,'SUBJ####','EEG','AMICA'). This script will 
%automatically find pre existing subfolders (expected increasing numbers) 
%and create a new %subfolder so one subject can have many AMICA 
%decompositions available. You just need to keep track of the unique number
%for each amica run
%
%avgRefPCAReduction is a number. It can be 1 or 2 or maybe even more.
%Example: set to 3 if using EEG, EMG, and Noise recordings all separately
%avg ref to each other. 
%
% emailStr is a string with your email if you want to be notified when the
% hipergator ran your data
%
%
% Modified by Chang - change paths and memory requirement(08-26-2021)

% Generate Float File for AMICA
jobName = fileNameNoExt; %jobName is a string that hipergator shows when you run squeue or get gemails

%% Parameters of interest
pcaOverride = 0;
max_iter = 2500; %normally 2000?
max_duration = hours(2); %minutes, days, etc.

num_models = 1;
%NOTE: share_comps might be on! It likes to make AMICA crash.

%% Get file, make float
subFolderExists = 1;
subDirNum = 0;
while subFolderExists == 1
    subDirNum = subDirNum + 1;
    subFolderExists = ~isempty(dir(fullfile(amicaOutputFolder_local,num2str(subDirNum))));
end
display(subDirNum);

AMICAdir_local = fullfile(amicaOutputFolder_local,num2str(subDirNum));

% AMICAdir_unix = ['/blue/dferris/share/MindInMotion/Data/',subj,'/EEG/AMICA/',num2str(subDirNum),'/'];
AMICAdir_unix = [amicaOutputFolder_unix,num2str(subDirNum),'/'];

num_chans = EEG.nbchan;
num_frames = length(EEG.times);

%% Save float file for running ICA

tmpdata = EEG.data;
disp('Converting data to double...');
tmpdata = double(tmpdata);

[SUCCESS,MESSAGE,MESSAGEID] = mkdir(AMICAdir_local);
%write file (make sure eeglab writes set and float not just set or you get error)
EEG = pop_saveset( EEG, 'filepath', AMICAdir_local, 'filename', fileNameNoExt);


%% Save AMICA parameter file 
hipergator_dir = AMICAdir_unix;
floatFileName = [fileNameNoExt,'.fdt'];
% num_nodes = 2; % How many nodes to request
% num_tasks = 8; % Number of MPI jobs
% num_procs = 4; % number of cores per job
num_mem = 600; % memory per cpu used (in MB)
num_nodes = 64; % How many nodes to request
num_tasks = num_nodes; % Number of MPI jobs
num_procs = 1; % number of cores per job
num_mem = 600*num_models; % memory per cpu used (in MB) %originally it was set to 256MB but the nodes tend to be out-of-memory
num_tasks_per_node = 1; %ryan edit: noticed the code ran best when you only used a single CPU from each node so i would generally set this to 1

max_duration.Format = 'hh:mm:ss';
num_time = char(max_duration); % Wall time hh:mm:ss
% num_time = '24:00:00'; % Wall time hh:mm:ss

qos_name = 'dferris-b'; % Use 'dferris' or 'dferris-b' (burst allocation)
if pcaOverride ~= 0
    pca_keep = min([pcaOverride, rank(tmpdata), EEG.nbchan-avgRefPCAReduction]);
else
    pca_keep = min([rank(tmpdata), EEG.nbchan-avgRefPCAReduction]); % PCA value to use
end


display(['Num ch = ', num2str(EEG.nbchan)]);
display(['PCA Reduction = ', num2str(pca_keep)]);
fid = fopen(fullfile(AMICAdir_local, 'input_hipergator.param'),'w');

fprintf(fid,['files ' hipergator_dir floatFileName '\n']);
fprintf(fid,['outdir ' hipergator_dir '\n']);
fprintf(fid,['num_models %d\n'],num_models); %ryan temp switched to 3 
fprintf(fid,['num_mix_comps 3\n']);
fprintf(fid,['pdftype 0\n']);
fprintf(fid,['block_size 128\n']);
fprintf(fid,['max_iter %d\n'],max_iter); 
fprintf(fid,['num_samples 1\n']);
fprintf(fid,['data_dim %d\n'],num_chans);
fprintf(fid,['field_dim %d\n'],num_frames);
fprintf(fid,['field_blocksize 1\n']);
fprintf(fid,['share_comps 0\n']); %ryan temp switched to 1 %SHARE COMPS HERE
fprintf(fid,['share_start 100\n']);
fprintf(fid,['comp_thresh 0.990000\n']);
fprintf(fid,['share_iter 100\n']);
fprintf(fid,['lrate 0.100000\n']);
fprintf(fid,['minlrate 1.000000e-08\n']);
fprintf(fid,['lratefact 0.500000\n']);
fprintf(fid,['rholrate 0.050000\n']);
fprintf(fid,['rho0 1.500000\n']);
fprintf(fid,['minrho 1.000000\n']);
fprintf(fid,['maxrho 2.000000\n']);
fprintf(fid,['rholratefact 0.500000\n']);
fprintf(fid,['kurt_start 3\n']);
fprintf(fid,['num_kurt 5\n']);
fprintf(fid,['kurt_int 1\n']);
fprintf(fid,['do_newton 1\n']);
fprintf(fid,['newt_start 50\n']);
fprintf(fid,['newt_ramp 10\n']);
fprintf(fid,['newtrate 1.000000\n']);
fprintf(fid,['do_reject 0\n']);
fprintf(fid,['numrej 15\n']);
fprintf(fid,['rejsig 3.000000\n']);
fprintf(fid,['rejstart 1\n']);
fprintf(fid,['rejint 1\n']);
%note: there can be memory issues if your dataset is large and num_tasks is
%greater than 1. You can either 1) use PCA reduction on the data or 2) set 
%max_threads to 1 and set num_tasks to be a large number with one task per 
%node. AMICA runs quickly this way
% fprintf(fid,['max_threads %d\n'],num_tasks); 
fprintf(fid,['max_threads %d\n'],1); % max_threads = 1 run locally
fprintf(fid,['writestep 250\n']); %ryan switched from 10 to 250 since performance was slowing down a lot during writing
fprintf(fid,['write_nd 0\n']);
fprintf(fid,['write_LLt 1\n']); %ryan switched from 0 to 1
fprintf(fid,['decwindow 1\n']);
fprintf(fid,['max_decs 3\n']);
fprintf(fid,['update_A 1\n']);
fprintf(fid,['update_c 1\n']);
fprintf(fid,['update_gm 1\n']);
fprintf(fid,['update_alpha 1\n']);
fprintf(fid,['update_mu 1\n']);
fprintf(fid,['update_beta 1\n']);
fprintf(fid,['invsigmax 100.000000\n']);
fprintf(fid,['invsigmin 0.000000\n']);
fprintf(fid,['do_rho 1\n']);
fprintf(fid,['load_rej 0\n']);
fprintf(fid,['load_W 0\n']);
fprintf(fid,['load_c 0\n']);
fprintf(fid,['load_gm 0\n']);
fprintf(fid,['load_alpha 0\n']);
fprintf(fid,['load_mu 0\n']);
fprintf(fid,['load_beta 0\n']);
fprintf(fid,['load_rho 0\n']);
fprintf(fid,['load_comp_list 0\n']);
fprintf(fid,['do_mean 1\n']);
fprintf(fid,['do_sphere 1\n']);
fprintf(fid,['doPCA 1\n']);
fprintf(fid,['pcakeep %d\n'],pca_keep);
fprintf(fid,['pcadb 30.000000\n']);
fprintf(fid,['byte_size 4\n']);
fprintf(fid,['doscaling 1\n']);
fprintf(fid,['scalestep 1\n']);
fclose(fid);

%% Save hipergator job script
fid = fopen(fullfile(AMICAdir_local, 'runAMICA_hipergator.sh'),'w');
fprintf(fid,['#!/bin/sh\n']);
fprintf(fid,['#SBATCH --job-name=%s_ICA # Job name\n'],jobName);
fprintf(fid,['#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n']);
fprintf(fid,['#SBATCH --mail-user=%s  # Where to send mail\n'],emailStr);	
fprintf(fid,['#SBATCH --nodes=%d                    # Use one node\n'],num_nodes);
fprintf(fid,['#SBATCH --ntasks=%d                   # Run a single task\n'],num_tasks);	
fprintf(fid,['#SBATCH --ntasks-per-node=%d          # number of tasks per node\n'],num_tasks_per_node); %ryan edit: added this to have more control over how tasks are split across nodes
fprintf(fid,['#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n'],num_procs);
fprintf(fid,['#SBATCH --mem-per-cpu=%dmb                  # Total memory limit\n'],num_mem);
fprintf(fid,['#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node\n']);
fprintf(fid,['#SBATCH --time=%s              # Time limit hrs:min:sec\n'],num_time);
fprintf(fid,['#SBATCH --output=%sAmicaOUT.out     # Standard output and error log\n'], AMICAdir_unix);
fprintf(fid,['#SBATCH --account=dferris	     # Account name\n']);
fprintf(fid,['#SBATCH --qos=%s		     # Quality of service name\n\n'],qos_name);
% fprintf(fid,['#SBATCH --constraint=haswell|skylake|opteron|dhabi|broadwell|sandy-bridge']); %added 9/6/2019 because some nodes werent up to date and would throw an error
fprintf(fid,['#SBATCH --partition=hpg2-compute']); %added 12/15/2020 because recommended by hipergator IT when random nodes on older part of cluster would cause AMICA to crash. This limits selection to hipergator 2
fprintf(fid,['# Run your program with correct path and command line options\n']);
fprintf(fid,['module load ufrc\n']);
% fprintf(fid,['module load intel openmpi\n']);
fprintf(fid,['module load intel/2020 openmpi/4.0.3\n']); %changed 2020-04-14 based on Amanda's suggestion
% fprintf(fid,['mpirun /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sParticipant%s_input_hipergator.param\n'],hipergator_dir, num2str(Participant));
% fprintf(fid,['srun --mpi=pmix /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir);
% fprintf(fid,['srun --mpi=pmix_v2 /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir);
% fprintf(fid,['srun --mpi=pmix_v2 /ufrc/dferris/share/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir);
%fprintf(fid,['srun --mpi=pmix_v3 /blue/dferris/share/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir); %switched to v3 on 2020-04-14 based on Amanda's suggestion
fprintf(fid,['srun --mpi=pmix_v3 /blue/dferris/share/s.peterson/test/AMICA_15/amica15ub ' hipergator_dir 'input_hipergator.param\n']); %switched to v3 on 2020-04-14 based on Amanda's suggestion

fclose(fid);

disp('Done creating AMICA stuff');
disp('Use MobaXTerm to submit your job (use right click to paste)');
disp(newline);
disp(['cd ' AMICAdir_unix]);
disp(['sbatch runAMICA_hipergator.sh']);
disp('squeue -A dferris');

txt = ['cd ' AMICAdir_unix, ';' 'sbatch runAMICA_hipergator.sh;'];
disp(txt)
end