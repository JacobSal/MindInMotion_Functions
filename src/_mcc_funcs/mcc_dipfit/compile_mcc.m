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
%%
% global ADD_DIPFIT_COMPILE_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_DIPFIT_COMPILE_SUBMODS = true;
%## Add Study, Src, & Script Paths
path(SRC_DIR,path);
path(STUDY_DIR,path);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%% (METHOD 1) HARDCODE DEPENDENCIES
set_workspace
% fprintf('Adding necessary paths to MATLAB path\n');
% PATHS = [];
% PATHS.submods_dir = [fileparts(SRC_DIR) filesep 'submodules'];
% PATHS.functions_dir = [SRC_DIR filesep '_functions'];
dep = {SCRIPT_DIR,...
    [PATHS.submods_dir filesep 'fieldtrip/contrib/spike'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/fileexchange'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/signal'],...    
    [PATHS.functions_dir filesep '_override/simbio'],... % [PATHS.submods_dir filesep 'fieldtrip/external/simbio'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/stats'],...
    [PATHS.submods_dir filesep 'fieldtrip/fileio'],...
    [PATHS.submods_dir filesep 'fieldtrip/forward'],...
    [PATHS.submods_dir filesep 'fieldtrip/inverse'],...
    [PATHS.submods_dir filesep 'fieldtrip/plotting'],...
    [PATHS.submods_dir filesep 'fieldtrip/preproc'],...
    [PATHS.submods_dir filesep 'fieldtrip/utilities'],...
    [PATHS.submods_dir filesep 'postAmicaUtility'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/freesurfer']};
%     [PATHS.submods_dir filesep 'eeglab']};
% %(12/05/2024) JS, removed "[PATHS.submods_dir filesep
% %'fieldtrip/external/simbio']" as I think its messing with parallel
% %processing
% cellfun(@(x) path(x,path),dep);
% %(12/05/2024) JS, doing path(x,path) instead of path(path,x) so that it
% %will prioritize these paths.
%##
%- these are files that you should include if the compiler has trouble
%finding dependencies simply based on pathing. 
files_to_compile={[SCRIPT_DIR filesep 'mcc_dipfit.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_sourcemodel.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_leadfield.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_dipolefitting.m'],...
    [PATHS.submods_dir filesep 'eeglab/plugins/dipfit/eeglab2fieldtrip.m'],...
    [PATHS.submods_dir filesep 'postAmicaUtility/pop_loadmodout.m'],...
    [PATHS.submods_dir filesep 'eeglab/functions/popfunc/pop_loadset.m'],...
    [PATHS.submods_dir filesep 'eeglab/functions/adminfunc/eeg_checkset.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/freesurfer/MRIread.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_headmodel.m']};
fprintf('Compiling...\n');
data_to_include={[SCRIPT_DIR filesep 'eeg_options.txt'],...
    [SCRIPT_DIR filesep 'eeg_optionsbackup.txt']};
%% (METHOD 2) AUTOMATICALLY DETERMINE DEPENDENCIES
%- these are files that you should include if the compiler has trouble
%finding dependencies simply based on pathing. 
%(12/11/2024) JS, this only works if the working path structure is setup to
%to ensure the right functions will be accessed during compiling (e.g., if
%you need to override a function this must be at the top of the path.)
%## METHOD 1
%- determine files
%{
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
ADD_DIPFIT_COMPILE_SUBMODS = true;
set_workspace
fout = matlab.codetools.requiredFilesAndProducts([SCRIPT_DIR filesep 'mcc_dipfit.m']);
par_save(fout,SCRIPT_DIR,'dependencies_openod.mat');
%}
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
%{
% ADD_DIPFIT_COMPILE_SUBMODS = true;
set_workspace
%- add files
files_to_compile = par_load(SCRIPT_DIR,'dependencies.mat');
keepinds = zeros(1,length(files_to_compile));
paths_to_add = cell(1,length(files_to_compile));
data_to_include = cell(1,length(files_to_compile));
% keepinds = [];
% paths_to_add = {}; %cell(length(files_to_compile),1);
% data_to_include = {}; %cell(length(files_to_compile),1);
extra = {};
for i = 1:length(files_to_compile)
    if ~ispc
        files_to_compile{i} = convertPath2UNIX(files_to_compile{i});
    else
        files_to_compile{i} = convertPath2Drive(files_to_compile{i});
    end
    tmp = strsplit(files_to_compile{i},filesep);
    if any(contains(tmp,'@')) || any(contains(tmp,'+')) || any(contains(tmp,'private'))
        fprintf('Skipping: %s\n',strjoin(tmp,filesep));
    else
        % paths_to_add = [paths_to_add, strjoin(tmp(1:end-1),filesep)];
        paths_to_add{i} = [strjoin(tmp(1:end-1),filesep)];
        regchk = regexp(files_to_compile{i},'(?<=\.)[^.]*$','match');
        if strcmp(regchk,'m')
            keepinds(i) = 1;
        elseif strcmp(regchk,'mat')
            % data_to_include = [data_to_include, files_to_compile(i)];
            data_to_include{i} = files_to_compile(i);
        elseif strcmp(regchk,'mexw64')
            if ~ispc
                tmp = strsplit(files_to_compile{i},'.');
                tmp{end}= 'mexa64';
                tmp = strjoin(tmp,'.');
                files_to_compile{i} = tmp;
            end
            disp(i);
            keepinds(i) = 1;
        end
    end
end
% files_to_compile = files_to_compile(logical(keepinds));
% files_to_compile = files_to_compile(~cellfun(@isempty,files_to_compile));
% data_to_include = data_to_include(~cellfun(@isempty,data_to_include));
paths_to_add = paths_to_add(~cellfun(@isempty,paths_to_add));
%## FINAL CHECK FOR DUPLICATES
files_to_compile = files_to_compile(logical(keepinds));
out = cellfun(@(x) strsplit(x,filesep),files_to_compile,'UniformOutput',false);
out_n = cellfun(@(x) regexp(x{end},'[^.]*(?=\.)','match'),out);
out_e = cellfun(@(x) regexp(x{end},'(?<=\.)[^.]*$','match'),out);
[~,ind,ic] = unique(out_n);
ind_keep = zeros(1,length(files_to_compile));
for o = 1:length(ind)
    tmp = out{ind(o)};
    t_n = out_n{ind(o)};
    t_e = out_e{ind(o)};
    ii = find(strcmp(t_n,out_n));
    if length(ii) > 1
        t_e = out_e(ii);
        i1 = strcmp(t_e,'m');
        if ~ispc
            i2 = strcmp(t_e,'mexa64');
        else
            i2 = strcmp(t_e,'mexw64');
        end
        if any(i2)
            % ind_keep = [ind_keep, ii(i2)];
            ind_keep(o) = ii(i2);
        else
            fprintf('Unable to reconcile multiple file matches for: %s\n',t_n);
            fprintf('Values: '); fprintf('%s, ',t_e{:}); fprintf('\n');
        end
    else
        % ind_keep = [ind_keep, ind(o)];
        ind_keep(o) = ind(o);
    end
end
%-
% ind_keep = ind_keep(ind_keep~=0);
%-
files_to_compile = files_to_compile(logical(ind_keep));
paths_to_add = [unique(paths_to_add)];
cellfun(@(x) path(x,path),paths_to_add);
%- sanity check
%{
files_to_compile = files_to_compile';
paths_to_add = paths_to_add';
data_to_include = data_to_include';
%}
%(12/05/2024) JS, doing path(x,path) instead of path(path,x) so that it
%will prioritize these paths.
%- files that need to be included for later recall (e.g., .mat and .txt files)
data_to_include=[data_to_include...
    [SCRIPT_DIR filesep 'eeg_options.txt'],...
    [SCRIPT_DIR filesep 'eeg_optionsbackup.txt']];
%}
%% COMPILE
%- cd to source directory
mkdir([SCRIPT_DIR filesep '_out']);
if length(data_to_include) > 1
    eval(['mcc -m ' strjoin(files_to_compile,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' strjoin(data_to_include,' -a ')...
        ' -R -singleCompThread -v'])
elseif length(data_to_include) == 1
    eval(['mcc -m ' strjoin(files_to_compile,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' data_to_include{1},...
        ' -R -singleCompThread -v'])
else
    eval(['mcc -m ' strjoin(files_to_compile,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'],...
        ' -R -singleCompThread -v'])
end
toc