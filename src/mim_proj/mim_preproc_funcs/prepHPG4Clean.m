function [SUBJSTRUCT,fid] = prepHPG4Clean(SUBJSTRUCT,scriptName,scriptPath,emailStr,dataFolder,varargin)
%prepHPG4Clean
%#ok<*NBRAK>
%       
%   IN:
%       ## REQUIRED
%       SUBJSTRUCT, STRUCT
%           is the strucutre containing pertainant parallel processing
%           information
%
%       scriptName, CHAR
%            
%       scriptPath, CHAR
%            
%       emailStr, CHAR
%           
%       dataFolder, CHAR
%           
%       ## OPTIONAL
%       analysisTitle, CHAR
%           
%       logFileDir, CHAR
%            
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
filesep = '/'; % UNIX overrideuserN = 'jsalminen';
userN = 'jsalminen';
analysisTitle = 'JS_test';
logFileDir = ['M:' filesep userN filesep 'GitHub' filesep 'connectivityMIM',...
             filesep 'src' filesep '_data' filesep '_hpgOutput' filesep scriptName];
logFileDir = convertPath2UNIX(logFileDir,'dferris');
Defaults = {logFileDir,...
            analysisTitle};
p = inputParser;
%## REQUIRED
addRequired(p,'SUBJSTRUCT',@isstruct)
addRequired(p,'scriptName',@ischar)
addRequired(p,'scriptPath',@ischar)
addRequired(p,'emailStr',@ischar)
addRequired(p,'dataFolder',@ischar)
%## OPTIONAL
addOptional(p,'analysisTitle',Defaults{2},@ischar)
addOptional(p,'logFileDir',Defaults{1},@ischar)
parse(p, SUBJSTRUCT, scriptName, scriptPath, emailStr, dataFolder, varargin{:});
%## SET DEFAULTS
logFileDir = p.Results.logFileDir;
analysisTitle = p.Results.analysisTitle;
%% --------------------------------------------------------------------- %%
%## create an extra subdirectory if one already exists
for i = 1:length(SUBJSTRUCT)
    SUBJSTRUCT(i).EEGstats.hpgCleanPath = [];
    SubjF = [dataFolder SUBJSTRUCT(i).SubjectStr filesep 'tmpEEG' filesep '_hgout' filesep];
    subFolderExists = 1;
    subDirNum = 0;
    while subFolderExists == 1
        subDirNum = subDirNum + 1;
        subFolderExists = ~isempty(dir(fullfile(SubjF, num2str(subDirNum))));
    end
    mkdir(fullfile(SubjF, num2str(subDirNum)));
    SUBJSTRUCT(i).EEGstats.hpgCleanPath = fullfile(SubjF, num2str(subDirNum));
end
%## HARDCODE PARAMS
numNodes = 1;
numTasks = 1;
% numTasksPerNode = 1;
numProcs = 16; % 32 cores cause why not.
numMem = [num2str((length(SUBJSTRUCT)*1024*5)/numProcs) 'mb']; % changed to mem-per-cpu (06/24/2022) 50?40 GB seems to work well for 5 subjects worth of MIM data
numTime = '04:00:00'; % starting with 4hrs, but may need to be changed
qosName = 'dferris';
moduleName = 'matlab/2020b';
%% RECONFIGURE PATH FOR UNIX
mkdir(scriptPath)
%% Save hipergator job script
%if exist delete
fid = fopen(fullfile(scriptPath, sprintf('run_%s.sh',scriptName)),'w');
fprintf(fid,['#!/bin/bash\n']);
fprintf(fid,['#SBATCH --job-name=%s_%s # Job name\n'],analysisTitle,scriptName);
fprintf(fid,['#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)\n']);
fprintf(fid,['#SBATCH --mail-user=%s # Where to send mail\n'],emailStr);	
fprintf(fid,['#SBATCH --nodes=%d # Use one node\n'],numNodes);
fprintf(fid,['#SBATCH --ntasks=%d # Run a single task\n'],numTasks);	
% fprintf(fid,['#SBATCH --ntasks-per-node=%d # number of tasks per node\n'],numTasksPerNode); %ryan edit: added this to have more control over how tasks are split across nodes
fprintf(fid,['#SBATCH --cpus-per-task=%d # Number of CPU cores per task\n'],numProcs);
fprintf(fid,['#SBATCH --mem-per-cpu=%s# Total memory limit\n'],numMem);
fprintf(fid,['#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node\n'],'');
fprintf(fid,['#SBATCH --time=%s # Time limit hrs:min:sec\n'],numTime);
% fprintf(fid,['## NOTE: %u=userID, %x=jobName, %N=nodeID, %j=jobID, %A=arrayID, %a=arrayTaskID\n']);
fprintf(fid,['#SBATCH --output=%s' filesep '%s_%s.out # Standard output\n'],logFileDir,analysisTitle,scriptName);
fprintf(fid,['#SBATCH --error=%s' filesep '%s_%s.err # error log\n'],logFileDir,analysisTitle,scriptName);
fprintf(fid,['#SBATCH --account=%s # Account name\n'],qosName);
fprintf(fid,['#SBATCH --qos=%s-b # Quality of service name\n'],qosName);
fprintf(fid,['#SBATCH --partition=hpg-default # cluster to run on, use slurm command ''sinfo -s''\n']); %added 12/15/2020 because recommended by hipergator IT when random nodes on older part of cluster would cause AMICA to crash. This limits selection to hipergator 2
fprintf(fid,['# Run your program with correct path and command line options\n']);
fprintf(fid,'\n\n');
fprintf(fid,['echo "Date              = $(date)"\n',...
             'echo "Hostname          = $(hostname -s)"\n',...
             'echo "Working Directory = $(pwd)"\n',...
             'echo ""\n',...
             'echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"\n',...
             'echo "Number of Tasks Allocated      = $SLURM_NTASKS"\n',...
             'echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"\n\n']);
fprintf(fid,('module load %s\n\n'),moduleName);
% fprintf(fid,('cd %s\n'),convertPath2UNIX(scriptPath,'dferris'));
% fprintf(fid,('matlab -nodisplay -r ''%s; exit;'''),scriptName);
fprintf(fid,['# Create a temporary directory on scratch\n']);
fprintf(fid,'mkdir -p ./$SLURM_JOB_ID\n\n');
fprintf(fid,'# Kick off matlab\n');
fprintf(fid,'matlab -nodisplay < %s\n\n',[convertPath2UNIX(scriptPath,'dferris') filesep scriptName '.m']);
fprintf(fid,'# Cleanup local work directory\n');
fprintf(fid,'rm -rf ./$SLURM_JOB_ID\n');
fclose(fid);

disp('Done creating SLURM script');
disp('Use MobaXTerm/CommandPrompt to submit your job by typing the ''ssh'' command below');
disp('ssh <GatorLink Username>@hpg.rc.ufl.edu');
disp(newline);
disp(['cd ' convertPath2UNIX(scriptPath,'dferris')]);
fprintf(['sbatch run_%s.sh\n'],scriptName);
disp('squeue -A dferris');

% commandout = ['sbatch ' AMICAdir_unix 'runAMICA_hipergator.sh'];
%##
toc
end