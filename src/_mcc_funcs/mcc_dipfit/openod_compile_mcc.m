%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 
%- compile exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/_compiled/mcc_dipfit/run_compile_mcc.sh
%- run exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/c_run_mim_mcc_dipfit.sh
%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(fileparts(SCRIPT_DIR)); % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = fileparts(fileparts(SCRIPT_DIR)); % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%% SET WORKSPACE
%## Add Study, Src, & Script Paths
path(SRC_DIR,path);
path(STUDY_DIR,path);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
ADD_DIPFIT_COMPILE_SUBMODS = true;
set_workspace
%% DETERMINE FILES?
% fout = matlab.codetools.requiredFilesAndProducts([SCRIPT_DIR filesep 'mcc_dipfit.m']);
% par_save(fout,SCRIPT_DIR,'dependencies_openod.mat');
%% COMPILE
%-
fout={[SCRIPT_DIR filesep 'mcc_dipfit.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_sourcemodel.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_leadfield.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_dipolefitting.m'],...
    [PATHS.submods_dir filesep 'eeglab/plugins/dipfit/eeglab2fieldtrip.m'],...
    [PATHS.submods_dir filesep 'postAmicaUtility/pop_loadmodout.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_headmodel.m']};
%- files that need to be included for later recall (e.g., .mat and .txt files)
data_to_include={[SCRIPT_DIR filesep 'eeg_options.txt'],...
    [SCRIPT_DIR filesep 'eeg_optionsbackup.txt']};
mkdir([SCRIPT_DIR filesep '_out']);
eval(['mcc -m ' strjoin(fout,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' strjoin(data_to_include,' -a ')...
        ' -R -singleCompThread -v'])
toc