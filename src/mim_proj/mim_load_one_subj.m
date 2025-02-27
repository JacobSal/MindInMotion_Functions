function [EEG] = mim_load_one_subj(fname,fpath,subj_char,cond_char,group_char,sess_char,varargin)
%MIM_CREATE_ALLEEG Summary of this function goes here
%   Detailed explanation goes here
%     MATLAB CODE DOCUMENTATION

%     Code Designer: Jacob Salminen
%     Date: 11/25/2022
% 
%     Functions:
%     1. fem_eeglab_dipfit: Configures dipole fitting parameters for Finite Element Model (FEM) dipfit analysis.
%     2. custom_update_chanlocs: Updates EEG channel locations using custom electrode locations.
%     3. bem_eeglab_dipfit: Fits `dipole models to EEG data using Boundary Element Model (BEM) dipfit analysis.
% 
%     See individual function descriptions for more details on their functionality.
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%- logo
cat_logo();
%- time
tt = tic;

%## PATHS
if ~ispc
    pp = strsplit(path,':');
else
    pp = strsplit(path,';');
end
%-
tmp = regexp(pp,'fieldtrip');
tmp = pp(~cellfun(@isempty,tmp));
b1 = regexp(tmp,'fieldtrip','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_fieldtrip = b2{1}(1:b1{1});
    fprintf('fieldtrip path: %s\n',path_fieldtrip);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('fieldtrip path not found.\n');
    end
end
%-
tmp = regexp(pp,'eeglab');
tmp = pp(~cellfun(@isempty,tmp));
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_eeglab = b2{1}(1:b1{1});
    fprintf('EEGLAB path: %s\n',path_eeglab);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('EEGLAB path not found.\n');
    end
end

%## TIME
tic

%## DEFINE DEFAULTS
eeglab_bem_fpath  = [path_eeglab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
DEF_LOAD_STRUCT = struct('do_bem_dipfit',false,...
    'do_load_dipfit',true,...
    'dipfit_fpath',{''},...
    'chanloc_fpath',{''},...
    'eeglab_bem_fpath',[path_eeglab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep],...
    'mni_mri_fpath',[eeglab_bem_fpath filesep 'standard_mri.mat'],...
    'mni_vol_fpath',[eeglab_bem_fpath filesep 'standard_vol.mat'],...
    'mni_chan_fpath',[eeglab_bem_fpath filesep 'elec' filesep 'standard_1005.elc'],...
    'mni_coord_transform',[0 0 0; 0 0 -1.5708; 1 1 1],...
    'dip_num',1,...
    'dip_plot','off', ...
    'recalc_iclabel',true);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'fname',@ischar);
addRequired(p,'fpath',@ischar);
addRequired(p,'subj_char',@ischar);
addRequired(p,'cond_char',@ischar);
addRequired(p,'group_char',@ischar);
addRequired(p,'sess_char',@ischar);
%## PARAMETER
addParameter(p,'LOAD_STRUCT',DEF_LOAD_STRUCT,@(x) validate_struct(x,DEF_LOAD_STRUCT));
parse(p,fname,fpath,subj_char,cond_char,group_char,sess_char,varargin{:});
%## SET DEFAULTS
LOAD_STRUCT = p.Results.LOAD_STRUCT;
%-
LOAD_STRUCT = set_defaults_struct(LOAD_STRUCT,DEF_LOAD_STRUCT);
%% ===================================================================== %%
%## CREATE STUDY
fprintf(1,'==== Creating Study ====\n')
fprintf(1,'Loading Subject %s\n',subj_char)
[~,EEG,~] = eeglab_loadICA(fname,fpath);
%## CHECKS
fprintf('%s) dipfit available: %i\n',EEG.subject,isfield(EEG,'dipfit'));
if LOAD_STRUCT.do_bem_dipfit
    fprintf('Trying Boundary Element Model Dipfit...\n');
    try
        EEG.dipfit.model;
        EEG.dipfit.coord_transform;
        EEG.dipfit.mrifile;
        EEG.dipfit.hdmfile;
        EEG.dipfit.coordformat;
    catch e
        fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
        EEG = bem_eeglab_dipfit(EEG,LOAD_STRUCT);
    end
end
%- eeg dipfit
if LOAD_STRUCT.do_load_dipfit
    fprintf('Trying to Load Dipfit...\n');
    out = par_load(LOAD_STRUCT.dipfit_fpath);
    EEG.dipfit = out;
    %- append important information
    [EEG] = fem_eeglab_dipfit(EEG,LOAD_STRUCT);
end
%- eeg channel locations
if ~isempty(LOAD_STRUCT.chanloc_fpath)
    fprintf('Attaching customelectrode locations...\n');
    try
        EEG = custom_update_chanlocs(EEG,LOAD_STRUCT);
    catch e
        fprintf('Error. couldn''t attach custom electrode configuration. Make sure your .mat file is formatted appropriately...\n%s\n',getReport(e));
    end
end
%## CHECK ICLABEL
if ~isfield(EEG.etc,'ic_classification') || LOAD_STRUCT.recalc_iclabel
    EEG = iclabel(EEG, 'lite');
end
%## MAKE EEG STRUCTURES FOR STUDY
EEG.subject = subj_char;
EEG.filepath = fpath;
EEG.filename = fname;
EEG.group = group_char;
EEG.condition = cond_char;
EEG.session = sess_char;
EEG = eeg_checkset(EEG,'eventconsistency');
fprintf('done. mim_load_one_subj.m: %0.2g s\n',toc(tt));
end
function [EEG] = fem_eeglab_dipfit(EEG,LOAD_STRUCT)
    % coord_trans = LOAD_STRUCT.mni_coord_transform;
    mni_mri = LOAD_STRUCT.mni_mri_fpath;
    mni_vol = LOAD_STRUCT.mni_vol_fpath;
    mni_chans = LOAD_STRUCT.mni_chan_fpath;
    tmp = [];
    if ~isfield(EEG.dipfit,'coord_transform')
        EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1]; %COORD_TRANSFORM_MNI;
        tmp = [tmp, 'added default coord_transform; '];
    end
    if ~isfield(EEG.dipfit,'mrifile')
        EEG.dipfit.mrifile = mni_mri;
        tmp = [tmp, 'added default mrifile; '];
    end
    if ~isfield(EEG.dipfit,'hdmfile')
        EEG.dipfit.hdmfile = mni_vol;
        tmp = [tmp, 'added default hdmfile; '];
    end
    if ~isfield(EEG.dipfit,'coordformat')
        EEG.dipfit.coordformat = 'MNI';
        tmp = [tmp, 'added default coordformat; '];
    end
    if ~isfield(EEG.dipfit,'chanfile')
        EEG.dipfit.chanfile = mni_chans;
        tmp = [tmp, 'added default chanfile; '];
    end
    if ~isfield(EEG.dipfit,'chansel')
        EEG.dipfit.chansel = (1:EEG.nbchan);
        tmp = [tmp, 'added default chansel; '];
    end
    EEG.dipfit.proc_details = tmp;
end
function [EEG] = custom_update_chanlocs(EEG,LOAD_STRUCT)
    mni_chans = LOAD_STRUCT.chanloc_fpath;
    tmp = load(mni_chans,'-mat');
    chanlocs_new = tmp.chanlocs_new;
    nodatchans_new = tmp.nodatchans_new;
    %- update the EEG electrode locations
    % Be cautious that not all electrodes are EEG
    % Sanity check: if we have 120 electrodes digitized
    fprintf('Found total of %i electrodes',length(chanlocs_new));
    for p = 1:length(chanlocs_new)
        elec_idx = find(strcmpi(chanlocs_new(p).labels,{EEG.chanlocs(:).labels}));
        if ~isempty(elec_idx)
            % update all available fields
            EEG.chanlocs(elec_idx).X = chanlocs_new(p).X;
            EEG.chanlocs(elec_idx).Y = chanlocs_new(p).Y;
            EEG.chanlocs(elec_idx).Z = chanlocs_new(p).Z;
            EEG.chanlocs(elec_idx).theta = chanlocs_new(p).theta;
            EEG.chanlocs(elec_idx).radius = chanlocs_new(p).radius;
            EEG.chanlocs(elec_idx).sph_theta = chanlocs_new(p).sph_theta;
            EEG.chanlocs(elec_idx).sph_phi = chanlocs_new(p).sph_phi;
        end
    end
    % Add fiducials location 
    if isempty(EEG.chaninfo.nodatchans)
        EEG.chaninfo.nodatchans = nodatchans_new;
    end
    EEG = eeg_checkchanlocs(EEG); % check the consistency of the chanloc structure
end

function [EEG] = bem_eeglab_dipfit(EEG,LOAD_STRUCT)
    coord_trans = LOAD_STRUCT.mni_coord_transform;
    mni_mri = LOAD_STRUCT.mni_mri_fpath;
    mni_vol = LOAD_STRUCT.mni_vol_fpath;
    mni_chans = LOAD_STRUCT.mni_chan_fpath;
    dip_num = LOAD_STRUCT.dip_num;
    do_dip_plot = LOAD_STRUCT.do_dip_plot;
    %## FIT DIPOLES TO HEADMODEL IF NEEDED
    fprintf('MNI pop_dipfit_settings...\n');
    %## Check Files
    tmp = [];
    if ~isfield(EEG.dipfit,'coord_transform')
        EEG.dipfit.coord_transform = coord_trans;
        tmp = [tmp, 'added default coord_transform; '];
    end
    if ~isfield(EEG.dipfit,'mrifile')
        EEG.dipfit.mrifile = mni_mri;
        tmp = [tmp, 'added default mrifile; '];
    end
    if ~isfield(EEG.dipfit,'hdmfile')
        EEG.dipfit.hdmfile = mni_vol;
        tmp = [tmp, 'added default hdmfile; '];
    end
    if ~isfield(EEG.dipfit,'coordformat')
        EEG.dipfit.coordformat = 'MNI';
        tmp = [tmp, 'added default coordformat; '];
    end
    if ~isfield(EEG.dipfit,'chanfile')
        EEG.dipfit.chanfile = mni_chans;
        tmp = [tmp, 'added default chanfile; '];
    end
    if ~isfield(EEG.dipfit,'chansel')
        EEG.dipfit.chansel = (1:EEG.nbchan);
        tmp = [tmp, 'added default chansel; '];
    end
    EEG.dipfit.comment = tmp;
   
    %- pop_multifit.m
    %- DIPFIT (see. ft_dipolefitting())
    fprintf('pop_multifit...\n');
    EEG = pop_multifit(EEG,[],'dipoles',dip_num,'dipplot',do_dip_plot);
%     dipfit = EEG.dipfit;
%     save([fPath filesep sprintf('%s_dipfit_mni.mat',EEG.subject)],'dipfit');
end