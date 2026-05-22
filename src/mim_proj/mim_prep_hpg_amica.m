function [EEG,cmd_out] = mim_prep_hpg_amica(EEG,fdt_fpath,amica_out_fpath,email_char,varargin)
%       EEG is the EEG struct you want to save on M drive and run AMICA on via the
%       hipergator
%
%   IN: 
%       EEG, STRUCT (EEGLAB)
%           is the EEG struct you want to save on M drive and run AMICA on via the
%           hipergator           
%       float_fPath, CHAR
%           is the path where your EEG float (.fdt) file is located. (hint:
%           its usually where your .set file is.
%       amica_out_fPath, CHAR
%           place where you would like your AMICA .params & .sh files to be
%           saved. (recommended: save them to the same folder as your
%           cleaned .set and .fdt file.
%       email_char, CHAR
%           is a string with your email if you want to be notified when the
%           hipergator ran your data
%
%   OUT: 
%       EEG, struct
%           
%       cmd_out, CHAR
%           a command line entry that runs the .sh file use srun & pmix_v3  
%## TIME
tic
%## DEFINE DEFAULTS
%-- precalculations for smoother automation
pcakeepcalc = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('%s) Num Chans: %i\nRank of Data:\nsum(eig(cov(double(EEG.data''))) > 1e-7) = %i\nsee. (Kim H, Luo J, Chu S, et al. (2023)) for calculation.\n', ...
    EEG.subject,EEG.nbchan,pcakeepcalc)
%(02/21/2025) JS, A paper on Ghost IC's may suggest setting "pcakeep" parameter
% in the ICA algorithm to this value for more robust decompositions:
% Kim, H., Luo, J., Chu, S., Cannard, C., Hoffmann, S., & Miyakoshi, M. 
% (2023). ICA’s bug: How ghost ICs emerge from effective rank deficiency caused by EEG electrode interpolation and incorrect re-referencing. 
% Frontiers in Signal Processing, 3. https://doi.org/10.3389/frsip.2023.1064138
DEF_AMICA_PARAMS = struct('num_models',1, ...
    'max_iter',2000, ... % (02/21/2025) JS, using recommendation from (Kim H, Luo J, Chu S, et al. (2023))
    'pcakeep',pcakeepcalc, ... % (02/21/2025) JS, using recommendation from (Kim H, Luo J, Chu S, et al. (2023))
    'num_chans',EEG.nbchan, ...
    'num_frames',length(EEG.times), ...
    'param_fname','input_hipergator.param', ...
    'minlrate',1e-6); % (02/21/2025) JS, using recommendation from (Kim H, Luo J, Chu S, et al. (2023))
DEF_BASH_PARAMS = struct('num_nodes',64, ...
    'num_tasks',64, ...
    'num_mem',ceil(512*1*1.5), ...
    'num_tasks_per_node',1, ...
    'max_duration','03:00:00', ...
    'qos_name','dferris-b', ...
    'account_name','dferris', ...
    'amica_dll_fpath',[filesep 'blue' filesep 'dferris' filesep, ...
        'share' filesep 's.peterson' filesep 'rhel9_amica_15' filesep,...
        'amica17ub'], ...
    'amica_bash_fname','run_amica_hipergator.sh');

p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'fdt_fpath',@ischar);
addRequired(p,'amica_out_fpath',@ischar);
addRequired(p,'email_char',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'AMICA_PARAMS',DEF_AMICA_PARAMS,@(x) validate_struct(x,DEF_AMICA_PARAMS));
addParameter(p,'BASH_PARAMS',DEF_BASH_PARAMS,@(x) validate_struct(x,DEF_BASH_PARAMS));
%-- parse
parse(p,EEG,fdt_fpath,amica_out_fpath,email_char,varargin{:});
%## SET DEFAULTS
AMICA_PARAMS = p.Results.AMICA_PARAMS;
AMICA_PARAMS = set_defaults_struct(AMICA_PARAMS,DEF_AMICA_PARAMS); 
BASH_PARAMS = p.Results.BASH_PARAMS;
BASH_PARAMS = set_defaults_struct(BASH_PARAMS,DEF_BASH_PARAMS); 
%(05/10/2026) JS, udpating to conform to a more efficient hyperocmputing setup.
AMICA_PARAMS.num_tasks = BASH_PARAMS.num_tasks;
%% ===================================================================== %%
%## STORE AMICA/BASH PARAMS
EEG.etc.amica_run = AMICA_PARAMS;
EEG.etc.amica_run.bash_params = BASH_PARAMS;
EEG.etc.amica_run.amica_out_fPath = amica_out_fpath;

%## CONVERT AMICA OUTPUT PATH
if ~exist(fdt_fpath,'file')
    error('ERROR. %s does not exist',fdt_fpath);
end
if ~ispc
    amica_out_fpath = convertPath2UNIX(amica_out_fpath);
else
    amica_out_fpath = convertPath2Drive(amica_out_fpath);
end
%-- make amica directory
if ~exist(amica_out_fpath,'file')
    fprintf('Making directory ''amica_out_fPath'' %s...\n',amica_out_fpath);
    mkdir(amica_out_fpath);
end

%## CREATE PARAM FILE
[~,paramf_fpath] = make_param_amica(fdt_fpath, amica_out_fpath, ...
    AMICA_PARAMS); %PCA_KEEP,...
    % EEG.nbchan, length(EEG.times), NUM_fMODELS, MAX_ITERS);

%## CREATE HIPERGATOR BASH FILE
[~,bash_out_fPath] = make_amica_bash(EEG.subject, amica_out_fpath, paramf_fpath, email_char, ...
    BASH_PARAMS);
% AMICA_DLL_PATH,...
% email_char, NUM_NODES, NUM_TASKS, NUM_TASKS_PER_NODE, NUM_MEM, NUM_TIME, QOS_NAME);
fprintf('Done creating .param and .sh files for AMICA.\n');
cmd_out = {['cd ' convertPath2UNIX(amica_out_fpath)];...
              sprintf('sbatch  %s',convertPath2UNIX(bash_out_fPath))};
              
%## SUBFUNCTION ======================================================== %%
function [fid,paramf_fpath] = make_param_amica(float_fPath, out_fPath, ...
        params)
    %## MAKE FILE
    paramf_fpath = [out_fPath filesep params.param_fname];
    fid = fopen(paramf_fpath,'w+'); %changed from 'w'
    %## ADD PARAMS
    fprintf(fid,'files %s\n',convertPath2UNIX(float_fPath));
    fprintf(fid,'outdir %s\n',convertPath2UNIX(out_fPath));
    fprintf(fid,'num_models %d\n',params.num_models); %ryan temp switched to 3 
    fprintf(fid,'num_mix_comps 3\n');
    fprintf(fid,'pdftype 0\n');
    fprintf(fid,'block_size 128\n');
    fprintf(fid,'max_iter %d\n',params.max_iter); 
    fprintf(fid,'num_samples 1\n');
    fprintf(fid,'data_dim %d\n',params.num_chans);
    fprintf(fid,'field_dim %d\n',params.num_frames);
    fprintf(fid,'field_blocksize 1\n');
    fprintf(fid,'share_comps 0\n'); 
    %(02/21/2025) RJD (JS), temp switched to 1
    %(02/21/2025) JS, not sure why he temporarily switched it to 1, but its
    %back to 0 as per a copied script.
    fprintf(fid,'share_start 100\n');
    fprintf(fid,'comp_thresh 0.990000\n');
    fprintf(fid,'share_iter 100\n');
    fprintf(fid,'lrate 0.100000\n');
    fprintf(fid,'minlrate %d\n',params.minlrate);
    %(02/21/2025) JS, A paper on Ghost IC's may suggest changing this minimum lambda rate to 
    %  a value greater than 10e-6. Value was originally 1.00000e-08:
    % Kim, H., Luo, J., Chu, S., Cannard, C., Hoffmann, S., & Miyakoshi, M. 
    % (2023). ICA’s bug: How ghost ICs emerge from effective rank deficiency caused by EEG electrode interpolation and incorrect re-referencing. 
    % Frontiers in Signal Processing, 3. https://doi.org/10.3389/frsip.2023.1064138
    fprintf(fid,'lratefact 0.500000\n');
    fprintf(fid,'rholrate 0.050000\n');
    fprintf(fid,'rho0 1.500000\n');
    fprintf(fid,'minrho 1.000000\n');
    fprintf(fid,'maxrho 2.000000\n');
    fprintf(fid,'rholratefact 0.500000\n');
    fprintf(fid,'kurt_start 3\n');
    fprintf(fid,'num_kurt 5\n');
    fprintf(fid,'kurt_int 1\n');
    fprintf(fid,'do_newton 1\n');
    fprintf(fid,'newt_start 50\n');
    fprintf(fid,'newt_ramp 10\n');
    fprintf(fid,'newtrate 1.000000\n');
    fprintf(fid,'do_reject 1\n'); 
    %(10/12/2021) RJD, turned on do_reject based on makoto processing pipeline recommendation
    % https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code#What_do_those_AMICA_parameters_mean.3F_.2804.2F06.2F2018_updated.29
    % "If you want to know which data points were rejected, you can check
    % EEG.etc.amica.Lht. If any datapoint shows 0, it means the datapoints were rejected by AMICA. 
    % Note also that you might want to use this info for further data cleaning if you don't mind rejecting randomly 
    % distributing single datapoints."
    fprintf(fid,'numrej 15\n');
    fprintf(fid,'rejsig 3.000000\n');
    fprintf(fid,'rejstart 1\n');
    fprintf(fid,'rejint 1\n');
    %(02/21/2025) RJD? (JS), note: there can be memory issues if your dataset 
    % is large and num_tasks is greater than 1. You can either 
    % 1) use PCA reduction on the data or 2) set 
    % max_threads to 1 and set num_tasks to be a large number with one task per 
    % node. AMICA runs quickly this way
    %(05/10/2026) JS, RJD was wrong with the revious comment and simply needed to 
    % use fewer nodes and tasks and increase number of threads. This value should be set
    % to 32 or 64
    % fprintf(fid,['max_threads %d\n'],num_tasks); 
    fprintf(fid,'max_threads %d\n',params.num_tasks); 
    fprintf(fid,'writestep 250\n'); 
    %(02/21/2025) RJD (JS), ryan switched from 10 to 250 since performance was slowing down a lot during writing
    fprintf(fid,'write_nd 0\n');
    fprintf(fid,'write_LLt 1\n'); 
    %(02/21/2025) RJD (JS), ryan switched from 0 to 1
    fprintf(fid,'decwindow 1\n');
    fprintf(fid,'max_decs 3\n');
    fprintf(fid,'update_A 1\n');
    fprintf(fid,'update_c 1\n');
    fprintf(fid,'update_gm 1\n');
    fprintf(fid,'update_alpha 1\n');
    fprintf(fid,'update_mu 1\n');
    fprintf(fid,'update_beta 1\n');
    fprintf(fid,'invsigmax 100.000000\n');
    fprintf(fid,'invsigmin 0.000000\n');
    fprintf(fid,'do_rho 1\n');
    fprintf(fid,'load_rej 0\n');
    fprintf(fid,'load_W 0\n');
    fprintf(fid,'load_c 0\n');
    fprintf(fid,'load_gm 0\n');
    fprintf(fid,'load_alpha 0\n');
    fprintf(fid,'load_mu 0\n');
    fprintf(fid,'load_beta 0\n');
    fprintf(fid,'load_rho 0\n');
    fprintf(fid,'load_comp_list 0\n');
    fprintf(fid,'do_mean 1\n');
    fprintf(fid,'do_sphere 1\n');
    fprintf(fid,'doPCA 1\n');
    fprintf(fid,'pcakeep %d\n',params.pcakeep);
    fprintf(fid,'pcadb 30.000000\n');
    fprintf(fid,'byte_size 4\n');
    fprintf(fid,'doscaling 1\n');
    fprintf(fid,'scalestep 1\n');
    fclose(fid);
end

%## SUBFUNCTION ======================================================== %%
function [fid,sh_out_fPath] = make_amica_bash(subj_char, out_fpath, param_fpath, email_char, ...
            params)
    sh_out_fPath = [out_fpath filesep params.amica_bash_fname];
    %## CREATE BASH FILE
    fid = fopen(sh_out_fPath,'w+');
    %## SHEBANG
    fprintf(fid,'#!/bin/bash\n'); 
    %(02/21/2025) JS, change from 'sh' to 'bash', supposedly better features?
    %## JOB NAME
    fprintf(fid,'#SBATCH --job-name=%s_AMICA\n',subj_char);
    %## MAIL AND STANDARD OUTPUT: Mail events (NONE, BEGIN, END, FAIL, ALL)
    fprintf(fid,'#SBATCH --mail-type=ALL\n');
    %-- Where to send mail
    fprintf(fid,'#SBATCH --mail-user=%s\n',email_char); 
    %-- Standard output and error log
    fprintf(fid,'#SBATCH --output=%s/%%j_amica_out.log\n', convertPath2UNIX(out_fpath)); 
    %## NODES, NTASKS, & MEMORY MGMT
    %-- Use one node
    fprintf(fid,'#SBATCH --nodes=%d\n',params.num_nodes);
    %-- Run a single task
    fprintf(fid,'#SBATCH --ntasks=%d\n',params.num_tasks);
    % number of tasks per node;
    fprintf(fid,'#SBATCH --ntasks-per-node=%d\n',params.num_tasks_per_node); 
    %(02/21/2025) RJD, added this to have more control over how tasks are split across nodes
    %-- Total memory limit
    fprintf(fid,'#SBATCH --mem-per-cpu=%dmb\n',params.num_mem);
    %-- Distribute tasks cyclically first among nodes and then among sockets within a node
    fprintf(fid,'#SBATCH --distribution=cyclic:cyclic\n'); 
    %## ALLOCATABLE TIME
    fprintf(fid,'#SBATCH --time=%s\n',params.max_duration);
    fprintf(fid,'#SBATCH --account=%s\n',params.account_name);
    fprintf(fid,'#SBATCH --qos=%s\n',params.qos_name);
    fprintf(fid,'#SBATCH --partition=hpg-default\n'); 
    %(12/15/2020) RJD? (JS), because recommended by hipergator IT when 
    % random nodes on older part of cluster would cause AMICA to crash. 
    % This limits selection to hipergator 2
    %## ECHO DETAILS & MODULE LOAD
    fprintf(fid,'cd %s\n',convertPath2UNIX(out_fpath));
    %--
    fprintf(fid,['echo "Date              = $(date)"\n',...
                 'echo "Hostname          = $(hostname -s)"\n',...
                 'echo "Working Directory = $(pwd)"\n',...
                 'echo ""\n',...
                 'echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"\n',...
                 'echo "Number of Tasks Allocated      = $SLURM_NTASKS"\n',...
                 'echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"\n',...
                 'echo "slurm_mem_per_cpu $SLURM_MEM_PER_CPU"\n',...
                 'echo "slurm_mem_per_gpu $SLURM_MEM_PER_GPU"\n',...
                 'echo "slurm_mem_per_node $SLURM_MEM_PER_NODE"\n']);
    %--
    fprintf(fid,'module purge\n'); 
    %(03/07/2023) JS, not sue if this is needed, but its use is encouraged on the HiperGator Wiki.
    fprintf(fid,'module load gcc/12.2.0 openmpi/5.0.7 openblas lapack mkl\n');
    fprintf(fid,'srun --mpi=pmix_v5 %s %s\n',convertPath2UNIX(params.amica_dll_fpath), convertPath2UNIX(param_fpath));
    fclose(fid);
end

function [struct_out] = set_defaults_struct(x,DEFAULT_STRUCT,varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'x')
addRequired(p,'DEFAULT_STRUCT');
%## OPTIONAL
addOptional(p,'print_chks',true,@islogical);
%##
parse(p,x,DEFAULT_STRUCT,varargin{:});
print_chks = p.Results.print_chks;
%% ===================================================================== %%
struct_out = x;
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
s_flags = zeros(length(fs2),1);
vals1 = struct2cell(x);
%-- check field value's class type
pchk_fun(sprintf('Setting structure defaults...\n'),print_chks)
for f = 1:length(fs2)
    if isstruct(DEFAULT_STRUCT.(fs2{f}))
        s_flags(f) = true;
    end
end
%-- 
for f = 1:length(fs2)
    ind1 = strcmp(fs2{f},fs1);
    if any(ind1)
        if isempty(vals1{ind1}) && ~ischar(vals1{ind1})
            pchk_fun(sprintf('Setting struct.%s to default value\n',fs2{f}),print_chks)
        
            struct_out.(fs1{ind1}) = DEFAULT_STRUCT.(fs2{f});
        elseif isstruct(DEFAULT_STRUCT.(fs2{f}))
            pchk_fun(sprintf('Recursive iteration on struct.%s...\n',fs2{f}),print_chks)
        
            tmp = set_defaults_struct(x.(fs2{f}),DEFAULT_STRUCT.(fs2{f}),print_chks);
            struct_out.(fs1{ind1}) = tmp;
        end
    else
        pchk_fun(sprintf('Setting struct.%s to default value\n',fs2{f}),print_chks)
        struct_out.(fs2{f}) = DEFAULT_STRUCT.(fs2{f});
    end
end
end

%##
function pchk_fun(str_in,print_chks)
    if print_chks
        fprintf(str_in);
    end
end


function [b] = validate_struct(x,DEFAULT_STRUCT,varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'x')
addRequired(p,'DEFAULT_STRUCT');
%## OPTIONAL
addOptional(p,'print_chks',true,@islogical);
%##
parse(p,x,DEFAULT_STRUCT,varargin{:});
print_chks = p.Results.print_chks;
%% ===================================================================== %%
b = false;
struct_name = inputname(2);
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
vals1 = struct2cell(x);
vals2 = struct2cell(DEFAULT_STRUCT);
%-- check field value's class type
pchk_fun(sprintf('Running structure checks...\n'),print_chks)
%-- check field names
chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
if ~all(chk)
    pchk_fun(sprintf('Remove field(s): %s\n',strjoin(fs1(~chk),',')),print_chks)
    pchk_fun(sprintf('Fields for struct do not match for %s\n',struct_name),print_chks)
    return
end
%- check field value's class type
for f = 1:length(fs2)
    ind = strcmp(fs2{f},fs1);
    if any(ind)
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            pchk_fun(sprintf('\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind})),print_chks)
            return
        end
    else
        pchk_fun(sprintf('Struct.%s is not present in input structure. Setting to default.\n',fs2{f}),print_chks)
    end
end
b = true;
end



end