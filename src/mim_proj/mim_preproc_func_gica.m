function [EEG,amica_cmd,CONFIG_STRUCT] = mim_preproc_func_gica(subj_char,save_fPath,fOutput,fDataset,varargin)
%MAIN_FUNC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%
% NOTES:
%       Code designed for Mind In Motion Study (NIHU01) ran from the
%       University of Florida.
%       prep step 1: High pass filter, cleanline, and merge
%       
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/23/2023, MATLAB 2019a
% Copyright (C) Chang Liu, Designer & Creator (20210820)
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu, Code Dweeb
% 
% Modified - 20230205 - add EMG high filter for ICC (Chang Liu)
% Modified - 20230223 - function wrapper (Jacob Salminen)

%## TIME
tic
%## DEFINE DEFAULTS
DEF_CONFIG_STRUCT = struct( ...
    'shareFolder', 'R:\Ferris-Lab\share\MindInMotion\Data\', ...
    'subjMgmtFolder', [fDataset filesep 'subject_mgmt' filesep 'Trial_Cropping_V2_test.xlsx'], ... %## FOLDER SETTINGS
    'fullRankAvRefBool', false, ... %## PREPROCESSING PARAMETERS
    'avgRefPCAReduction', 1, ...
    'std_threshold', 3.5, ...
    'CCA', struct( ... %## CCA PARAMETERS
        'lagAmount_samples', 1, ...
        'CCA_Rsq_thres', 0.2 ...
    ), ...
    'ASR_correct', struct( ... %## ASR PARAMETERS
        'burst_crit', 2, ...
        'arg_window', 0.3, ...
        'arg_flatline',1, ...
        'arg_hp',[0.25,0.95], ...
        'line_noise',3, ...
        'channel_crit',0.7 ...
    ), ...
    'chan_crit1', 0.70, ... %## CleanLine / channel rejection
    'wind_crit', 0.4, ...
    'chan_max_broken_time', 0.6, ...
    'chan_detected_fraction_threshold', 0.6, ...
    'flatline_crit', 'off', ...
    'line_noise_crit', 'off', ...
    'kurt_crit', 5, ...
    'chan_detect_num_iter', 5, ...
    'window_crit_tolerances', [-Inf, 10], ...
    'kurt_prob_std_rm', 0, ...
    'filter_lowCutoffFreqAMICA', 1.75, ... %## AMICA FILTERING
    'filter_AMICA_highPassOrder', 1650, ...
    'filter_highCutoffFreqAMICA', [], ...
    'filter_AMICA_lowPassOrder', [], ...
    'emailStr', 'jsalminen@ufl.edu' ... %## AMICA
);
% fullrankavgbool = default false for regular avg ref, true for full rank version
% avg_pca_red = 1 for just eeg, 3 for eeg+noise+emg % what is the difference?

%##
DEF_ICC_EEG_STRUCT = struct( ...
    'statsWindow', 4, ...
    'cleanWindow', [], ...
    'rhoSqThres_source', 0.65, ...
    'cleanXwith', 'X', ...
    'cleanYwith', 'no', ...
    'rerefX', 'no', ...
    'rerefY', 'no', ...
    'useBoundary', false, ...
    'useExternCalibData', false, ...
    'RTBool', false, ...
    'calcCCAonWholeData', false, ...
    'plotStatsOn', 1, ...
    'giveCleaningUpdates', 1, ...
    'plotVAF', 0, ...
    'filtXtype', 'no', ...
    'filtYtype', 'no' ...
);
% 2023-1-10 Use 0.65 for 4sec window %Default 0.9, NH3066 try 0.8, H3063 = 0.95, NH3068, H3072 try 0.85
% cleanWindow = empty here will cause the cleaning window length to match the stats window length exactly.
% statsWindow = length of stats window (that CCA is calculated on, in seconds)

%##
DEF_ICC_EMG_STRUCT = DEF_ICC_EEG_STRUCT;
DEF_ICC_EMG_STRUCT.rhoSqThres_source = 0.4;
%* (07/11/2023) JS, trying rhoSqThres_source 0.65
%* (07/12/2023) JS, trying rhoSqThres_source 0.9 as a less aggressive cutoff 
%* (10/30/2023) JS, trying rhoSqThres_source 0.6 as a  more aggressive cutoff
% (NA) Anon, with 4sec window, 0.9 is too high, 0 bad sources get removed...- -(EMG without highpass filter)

%##
do_cleanline = true;                           
do_processEMG = true;
%## Script Specific Params
cushion_sec = 5; %sec
%## PREPROCESSING PARAMS
do_EMG_highpass = 1; %use .set with highpass filter EMG
do_iCC = true;
%(02/10/2026) JS, trying to not use this for a run.
do_iCC_muscle = true;
do_CCA = false;
do_ASR = true;
do_autoreject_chans = true;
do_autoreject_chans_interp = false;
%(02/10/2026) JS, not using this for a run, but adding in a final average
%rereferncing post merge of all subjects.
cleaningMethod = '';
finalChToKeepForICA = {'EEG'};

%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'subj_char',@ischar);
addRequired(p,'fPath',@ischar);
addRequired(p,'fOutput',@ischar);
addRequired(p,'fDataset',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'CONFIG_STRUCT',DEF_CONFIG_STRUCT,@(x) validate_struct(x,DEF_CONFIG_STRUCT));
addParameter(p,'ICC_EEG_STRUCT',DEF_ICC_EEG_STRUCT,@(x) validate_struct(x,DEF_ICC_EEG_STRUCT));
addParameter(p,'ICC_EMG_STRUCT',DEF_ICC_EMG_STRUCT,@(x) validate_struct(x,DEF_ICC_EMG_STRUCT));

parse(p,subj_char,save_fPath,fOutput,fDataset,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETERS
CONFIG_STRUCT = p.Results.CONFIG_STRUCT;
CONFIG_STRUCT = set_defaults_struct(CONFIG_STRUCT,DEF_CONFIG_STRUCT);
ICC_EEG_STRUCT = p.Results.ICC_EEG_STRUCT;
ICC_EEG_STRUCT = set_defaults_struct(ICC_EEG_STRUCT,DEF_ICC_EEG_STRUCT);
ICC_EMG_STRUCT = p.Results.ICC_EEG_STRUCT;
ICC_EMG_STRUCT = set_defaults_struct(ICC_EMG_STRUCT,DEF_ICC_EMG_STRUCT);
%- MIMCONFIG
CONFIG_STRUCT.CleanEEG_output     = [fOutput filesep 'clean'];
CONFIG_STRUCT.amicaEEG_output     = [fOutput filesep 'clean'];
%--
autoChRejMethod = [num2str(CONFIG_STRUCT.std_threshold),'std'];
%% ===================================================================== %%
%## (SCRIPT PART 1)
EEG = [];
ALLEEG = [];
%## Find All Trials of Interest
fileList_TM = dir([save_fPath filesep 'TM*.set']);
fileList_SP = dir([save_fPath filesep 'SP*.set']);
fileList_Rest = dir([save_fPath filesep '*.set']);
inds = cellfun(@(x) contains(x,'rest','IgnoreCase',true),{fileList_Rest.name});
fileList_Rest = fileList_Rest(inds);
%(11/21/2024) JS, not including MotorImagery because not enough for all
%subjs
% fileList_Imagined = dir([save_fPath filesep 'MotorImagery*.set']);
% fileList = [fileList_Rest; fileList_TM; fileList_SP];% for this analysis, only include Rest, TM, and SP
fileList = fileList_Rest;
% if length(fileList) < 10
%     error('Subject %s missing condition',subj_char)
% end
%(02/19/2025) JS, removing this because of the new implementation of a
%linear mixed effects model which can handl subjects with partial
%conditions
%% Trial Level Processing
% The processing includess a highpass filter cutoff at 1Hz, and then use
% cleanline at 60Hz for all channels
% eeglab; %initialize EEG and ALLEEG as empty
for trial_i = 1:size(fileList,1)
    %Load trial
    EEG = pop_loadset('filename',fileList(trial_i).name,'filepath',save_fPath);
    %(11/22/2024) JS, trying to fix a bug where the IMU channels don't get
    %loaded if I load using the pop_loadset function. Actually they get
    %loaded but they are full of nan's not sure where the entire trace gets
    %removed
    %(02/21/2024) JS, fixed this, just an issue with some of the data
    %points, not all
    % tmp = load('-mat', [save_fPath filesep fileList(trial_i).name]);
    % tmp = tmp.EEG;
    % fid_eeg = fopen([save_fPath filesep tmp.data], 'r', 'ieee-le');
    % tmp.data    = fread(fid_eeg, [tmp.pnts*tmp.trials tmp.nbchan], 'float32')';
    % data = fread(fid_eeg, [tmp.nbchan tmp.trials*tmp.pnts], 'float32');
    % tmp.data = data;
    %##
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); %append to new EEG to ALLEEG
    %% Deal with trial cropping
    fprintf('Cropping trial %i for subject %s',trial_i,subj_char)
%     [DoCrop, ExactCrop, DoCrop_loadsol, ExactCrop_loadsol ] = CropTrialCheckFunc_checkLoadsol(subj_char,fileList(trial_i).name(1:end-4),...
%         mim_config.subjMgmtFolder);
    [ DoCrop, ExactCrop, DoCrop_loadsol, ExactCrop_loadsol ] = mim_check_trials(subj_char,fileList(trial_i).name(1:end-4),...
        CONFIG_STRUCT.subjMgmtFolder);
    if DoCrop_loadsol % 
        % remove gait events not in the ExactCrop_loadsol range
        % first convert the time --> latency
%             temp = [EEG.event.latency];
        disp('You are removing bad gait events');
        cushion_latency = EEG.srate*cushion_sec;
        crop_latency_range = EEG.srate*ExactCrop_loadsol + cushion_latency + 1;
        total_gaitevent_idx = find(strcmp('LHS',{EEG.event.type}) | strcmp('LTO',{EEG.event.type})|...
            strcmp('RHS',{EEG.event.type}) | strcmp('RTO',{EEG.event.type}));
        gaitevent_idx_inrange = [];
        for p = 1:size(crop_latency_range,1)
            gaitevent_idx_inrange = horzcat(gaitevent_idx_inrange,find( [EEG.event.latency] > ...
                crop_latency_range(p,1) & [EEG.event.latency] < crop_latency_range(p,2)));
        end
        [gaitevent_idx_outrange] = setdiff(total_gaitevent_idx,gaitevent_idx_inrange);
        EEG.event(gaitevent_idx_outrange) = [];

%             figure();plot(temp(gaitevent_idx_outrange)./500,ones(length(gaitevent_idx_outrange)) ,'o');hold on;plot(temp./500,zeros(length(temp)),'o');
%             figure();plot(temp(gaitevent_idx_inrange)./500,ones(length(gaitevent_idx_inrange)) ,'o');
    end
    if DoCrop %only deal with shortened trials for now (cannot ged rid of other events that happen during trial)
        disp('You are cropping the trial');
        startInd = find(strcmpi('TrialStart',{EEG.event.type}));
        startTime = EEG.times( round(EEG.event(startInd).latency) )/1000;
        endInd = find(strcmpi('TrialEnd',{EEG.event.type}));
        ExactCropLatencies = EEG.srate*ExactCrop(end)+EEG.event(startInd).latency;        
        endTime = EEG.times( round(ExactCropLatencies))/1000;
        %- add trial end
        EEG.event(end+1) = EEG.event(endInd);
        EEG.event(end).latency = ExactCropLatencies;
        EEG.event(end).datetime = [];

        cushion_sec = 5;
        ExactCrop_update = ExactCrop+cushion_sec;
        ExactCrop_update(1) = 0;
        ExactCrop_update(end,end) = endTime+cushion_sec;
        EEG = pop_select( EEG, 'time',ExactCrop_update); % This is wrong in the original code from Roehl. That code doesn't take into account the cushion time added for each trial
    end

    %## HP filter EEG, EMG, and Noise
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'chantype',{'EEG' 'EMG' 'Noise'}); %1 Hz HP filter (results in -6bB cutoff at 0.5 Hz)
    if do_processEMG
        EEG = pop_eegfiltnew(EEG, 'locutoff',20,'chantype',{ 'EMG' }); %high pass filter 20 Hz
    end
    %- Keep only a subset of chans
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    % back_imu_chans = find(contains({EEG.chanlocs.labels},'Back','IgnoreCase',true));
    % ls_chans = find(contains({EEG.chanlocs.labels},'GRF','IgnoreCase',true));
    %-- select subset of channels
    EEG = pop_select( EEG, 'channel', sort([eeg_chans,emg_chans,noise_chans])); 
    %-- select just EEG channels
    % EEG = pop_select( EEG, 'channel', sort([eeg_chans])); %
    %(02/10/2026) Don't need noise channels when not doing the gait stuff.    
    %-- redefine channels (numbering changed since we took subset of channels)
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); 
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    % back_imu_chans = find(contains({EEG.chanlocs.labels},'Back','IgnoreCase',true));
    % ls_chans = find(contains({EEG.chanlocs.labels},'GRF','IgnoreCase',true));
    %--
    %## Cleanline
    if do_cleanline % removal of line noise (set to remove 60 Hz below)
        % Reference before cleanline filter
        EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
        EEG = pop_cleanline(EEG, 'bandwidth',2,...
            'chanlist',[eeg_chans emg_chans noise_chans],...
            'computepower',1,...
            'linefreqs',[60 120],...
            'newversion',0,...
            'normSpectrum',0,...
            'p',0.01,...
            'pad',2, ...
            'plotfigures',0, ...
            'scanforlines',1, ...
            'sigtype','Channels', ...
            'taperbandwidth',2, ...
            'tau',100, ...
            'verb',1, ...
            'winsize',4, ...
            'winstep',1);
        [ALLEEG,~,~] = eeg_store(ALLEEG, EEG, CURRENTSET); %These params followed PREP pipeline
    else
        %Update ALLEEG with whatever changes you made to EEG variable before moving on to next trial
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET); %overwrites current set    
    end
end
%## Session level Processing
EEG = pop_mergeset(ALLEEG,1:length(ALLEEG), 0);
EEG.setname = sprintf('%s_merged_eeg',subj_char);
EEG.filename = sprintf('%s_merged_raw',subj_char);
[~, EEG, ~] = eeg_store( ALLEEG, EEG, 0 ); %append to new EEG to ALLEEG
%##
%## REMOVE UNNECESSARY ALLEEG FOR MEMORY
% ALLEEG = struct.empty;
%% ===================================================================== %%
%## SCRIPT PART 2
tic
fprintf('==== TIME REJECTIONS, ICANCLEAN, OR OTHER %s PREPROCESSING ====\n\n',subj_char) 
if ischar(finalChToKeepForICA) %e.g. just 'EEG', not {'EEG'} or {'EEG','Noise'}
    finalChToKeepForICA = {finalChToKeepForICA}; %turn 'EEG' into {'EEG'} to make code below more robust
end
%--
if CONFIG_STRUCT.fullRankAvRefBool
    avg_ref_pca_reduction = 0;
else
    avg_ref_pca_reduction = length(finalChToKeepForICA); %1 for just {'EEG'} and 3 for {'EEG','Noise','EMG'}
end

%% Load subject filtered merged file
if do_EMG_highpass
    fileName = [subj_char,'_HP1hz_cleanline_merge_EMG.set']; % edit here for file name
    cleaningMethod = horzcat(cleaningMethod,'EMG_HP');
else
    fileName = [subj_char,'_HP1hz_cleanline_merge.set']; % edit here for file name
end
% EEG = pop_loadset('filepath',mim_config.EEG_input,'filename',fileName);
EEG.urchanlocs = EEG.chanlocs; % keep one copy of old channel info; needed for chan_rej
%% Session Level Processing
EEG.setname = sprintf('%s_merged_eeg',subj_char);
EEG.filename = sprintf('%s_merged_raw',subj_char);
% eeglab redraw;
EEG.subject = subj_char;
%## KURTOSIS REJECTION OF CHANNELS AND REREFERENCING
%-- Re-ref EEG, EMG, and Noise to themselves
fprintf('==== %s 1.) REREFERENCING ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
%-- Reject bad channels
% pop_eegplot( EEG, 1, 1, 1);
if do_autoreject_chans
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    bad_eeg_chans = [];
    EEG.etc.chan_rej = [];
    fprintf('==== %s 2.) AUTOREJECTION OF BAD CHANS ====\n',subj_char)
    % testEEG = EEG;
    [EEG, badEEGch, badEMGch, badNoiseCh] = autoRejCh_func_CL(EEG,CONFIG_STRUCT.std_threshold,eeg_chans,emg_chans,noise_chans);

    bad_eeg_chans = cat(2,bad_eeg_chans,badEEGch);
    %-- Re-ref again (since remnants of rejected channels still exist    
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    %-- last average referencing
    fprintf('==== %s 3.) REREFERENCING pt.2 ====\n',subj_char)
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    %-- reject and re-ref again
    fprintf('==== %s 4.) AUTOREJECTION OF BAD CHANS pt.2 ====\n',subj_char)
    [EEG, badEEGch, badEMGch, badNoiseCh] = autoRejCh_func_CL(EEG,CONFIG_STRUCT.std_threshold,eeg_chans,emg_chans,noise_chans);

    bad_eeg_chans = cat(2,bad_eeg_chans,badEEGch);
    %-- Re-ref again (since remnants of rejected channels still exist 
    fprintf('==== %s 5.) REREFERENCING pt.3 ====\n',subj_char)
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    tt = strsplit(num2str(autoChRejMethod),'.');
    cleaningMethod = horzcat(cleaningMethod,[tt{1},'p',tt{2}]);
    EEG.etc.chan_rej = bad_eeg_chans;
end

%##
if do_autoreject_chans_interp
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    bad_eeg_chans = [];
    EEG.etc.chan_rej = [];
    tmp_eeg = EEG;
    fprintf('==== %s 2.) AUTOREJECTION OF BAD CHANS ====\n',subj_char)
    % testEEG = EEG;
    [tmp, badEEGch, badEMGch, badNoiseCh] = autoRejCh_func_CL(EEG,CONFIG_STRUCT.std_threshold,eeg_chans,emg_chans,noise_chans);
    EEG = tmp;
    bad_eeg_chans = cat(2,bad_eeg_chans,badEEGch);
    %-- Re-ref again (since remnants of rejected channels still exist from    
    % EEG = eeg_interp(EEG,badEEGch,'spherical',[EEG.times(1),EEG.times(end)]);
    % EEG = eeg_interp(EEG,badEEGch,'sphericalKang',[EEG.times(1),EEG.times(end)],[1e-8,3,50]);
    % EEG.data(badEEGch,:) = zeros(length(badEEGch),size(EEG.data,2));
    eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
    emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    %-- last average referencing
    fprintf('==== %s 3.) REREFERENCING pt.2 ====\n',subj_char)
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    %-- reject and re-ref again
    fprintf('==== %s 4.) AUTOREJECTION OF BAD CHANS pt.2 ====\n',subj_char)
    [tmp, badEEGch, badEMGch, badNoiseCh] = autoRejCh_func_CL(EEG,CONFIG_STRUCT.std_threshold,eeg_chans,emg_chans,noise_chans);
    EEG = tmp;
    bad_eeg_chans = cat(2,bad_eeg_chans,badEEGch);
    %-- Re-ref again (since remnants of rejected channels still exist from
    % EEG = eeg_interp(EEG,badEEGch,'sphericalKang',[EEG.times(1),EEG.times(end)],[1e-8,3,50]);
    % EEG = eeg_interp(EEG,badEEGch,'spherical',[EEG.times(1),EEG.times(end)]);
    %-- re-ref
    fprintf('==== %s 5.) REREFERENCING pt.3 ====\n',subj_char)
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    cleaningMethod = horzcat(cleaningMethod,autoChRejMethod);    
    %-- update channel types b/c that is always a good idea
    EEG = eeg_interp(tmp_eeg,bad_eeg_chans,'spherical',[tmp_eeg.times(1),tmp_eeg.times(end)]);
    cleaningMethod = horzcat(cleaningMethod,autoChRejMethod,'_interpc');
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    EEG.etc.chan_rej = bad_eeg_chans;
end
%% CLEANING ============================================================ %%

%## ICanClean Base
%(11/21/2024) JS, not using for find analysis
% pop_eegplot(EEG,1,0,0)
if do_iCC
    rho = num2str(ICC_EEG_STRUCT.rhoSqThres_source);
    % idx_rho = strfind(rho,'.');
    tt = strsplit(num2str(ICC_EMG_STRUCT.rhoSqThres_source),'.');
    cleaningMethod = horzcat(cleaningMethod,'_iCC',tt{1},'p',tt{2});
    fprintf('==== %s STARTING: ICANCLEAN ====\n',subj_char);
    [eeg_chans, emg_chans, noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,eeg_chans, noise_chans, 0, ICC_EEG_STRUCT);
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
end

%##
%(11/21/2024) JS, not using for find analysis
if do_iCC_muscle
    % rho_iCC_muscle = num2str(ICC_EMG_STRUCT.rhoSqThres_source);
    tt = strsplit(num2str(ICC_EMG_STRUCT.rhoSqThres_source),'.');
    cleaningMethod = horzcat(cleaningMethod,'_iCCEMG',tt{1}, 'p',tt{2});
    fprintf('==== %s STARTING: ICANCLEAN MUSCLE ====\n',subj_char);
    [eeg_chans, emg_chans, noise_chans] = getChannelTypes_func(EEG);
    [EEG] = iCanClean(EEG,[eeg_chans], [emg_chans], 0, ICC_EMG_STRUCT);
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
end

%## CCA
%(11/21/2024) JS, not using for find analysis
if do_CCA
    rho_CCA = num2str(CONFIG_STRUCT.CCA.CCA_Rsq_thres);
    idx_rho_CCA = strfind(rho_CCA,'.');
    cleaningMethod = horzcat(cleaningMethod,'_CCA','0p',rho_CCA(idx_rho_CCA+1:end));
    fprintf('==== %s STARTING: CCA ====\n',subj_char);
    [eeg_chans, emg_chans, noise_chans] = getChannelTypes_func(EEG)        
    EEG = autoLagCCA_wrap(EEG,CONFIG_STRUCT.CCA);
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
end

%## ASR
%(11/21/2024) JS, not using for find analysis
if do_ASR
    tmp_asr_params = CONFIG_STRUCT.ASR_correct;
    % [EEG_temp_clean,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,CONFIG_STRUCT); 
    %
    % BurstCriteria = num2str(CONFIG_STRUCT.ASR_correct.burstCritASR);
    % cleaningMethod = horzcat(cleaningMethod,'_ASRcorr',BurstCriteria);
    fprintf('==== %s STARTING: AUTOMATIC SUBSPACE RECONSTRUCTION ====\n',subj_char);
    % [EEG] = ASR_correction(EEG,CONFIG_STRUCT.ASR_correct);
    %## ASR w/ no reference values
    [eeg_chans, ~, ~] = getChannelTypes_func(EEG) ;
    EEG_only = pop_select(EEG,'channel',eeg_chans);
    % EEG_ASR = clean_artifacts(EEG_only, ...
    %     'FlatlineCriterion','off', ...
    %     'ChannelCriterion','off', ...
    %     'LineNoiseCriterion','off',...
    %     'Highpass','off', ...
    %     'BurstCriterion',tmp_asr_params.burstCritASR, ...
    %     'WindowCriterion','off',...
    %     'BurstRejection','off', ...
    %     'Distance','Euclidian');
    EEG_ASR = clean_rawdata(EEG_only, ...
        tmp_asr_params.arg_flatline, ...
        tmp_asr_params.arg_hp, ...
        tmp_asr_params.channel_crit, ...
        tmp_asr_params.line_noise, ...
        tmp_asr_params.burst_crit, ...
        tmp_asr_params.arg_window);
    %--
    % vis_artifacts(EEG_ASR,EEG_only);
    % figi = get(groot,'CurrentFigure');
    % figi.Position = [100,200,720,480];
    % exportgraphics(figi,[CONFIG_STRUCT.CleanEEG_output filesep 'rASR_out.png']);
    % close(figi);
    %--
    % fid = fopen([CONFIG_STRUCT.CleanEEG_output filesep 'info.txt'],'w');
    % fprintf(fid,'\n %.2f percent of frames were rejected\n', p_frames_rej);
    % % fprintf(fid,'\n %.2f channels were rejected\n', p_chan_rej);
    % fprintf(fid,'\n %.2f channels were rejected\n', tmp_eeg.nbchan-EEG.nbchan);
    % fclose(fid);
    %--
    % inds = cellfun(@(x) any(strcmp(x,{EEG_ASR.chanlocs.labels})),{EEG.chanlocs.labels});
    % EEG.data(inds,1:size(EEG_ASR.data,2)) = EEG_ASR.data;
    EEG = EEG_ASR;
    EEG_only = [];    
    EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
    %--
    wsv = strsplit(num2str(tmp_asr_params.arg_window),'.');
    csv = strsplit(num2str(tmp_asr_params.channel_crit),'.');
    cleaningMethod = horzcat(cleaningMethod,'_ASR', ...
        num2str(tmp_asr_params.burst_crit), ...
        ['w',wsv{1},'p',wsv{2}], ...
        ['c',csv{1},'p',csv{2}]);
end

%% Save EEG and log
EEG = rmfield(EEG,'icaweights');
EEG.icaweights = [];
EEG = rmfield(EEG,'icawinv');
EEG.icawinv = [];
EEG = rmfield(EEG,'icaact');
EEG.icaact = [];
EEG = rmfield(EEG,'icasphere');
EEG.icasphere = [];
EEG = rmfield(EEG,'icachansind');
EEG.icachansind = []; %also remove icachanind?
%--
preprocess_pipeline = [cleaningMethod];
fprintf('Saving %s EEG to %s\n',subj_char,fullfile(CONFIG_STRUCT.CleanEEG_output,preprocess_pipeline));
EEG.etc.CleanType  = preprocess_pipeline;
EEG.etc.Params     = CONFIG_STRUCT;
mkdir(CONFIG_STRUCT.CleanEEG_output)

% [eeg_chans,emg_chans,noise_chans] = getChannelTypes_func(EEG) ;
eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
biom_chans = find(strcmpi('BioM',{EEG.chanlocs.type}));
%(11/22/2024) JS, this may be a more robust way of grabbing channels, but
%I'm not going to use it as I perform this at the beginning.
% imu_chans = find(contains({EEG.chanlocs.labels},'Back','IgnoreCase',true));
% ls_chans = find(contains({EEG.chanlocs.labels},'GRF','IgnoreCase',true));
% re-ref
fprintf('==== %s FINAL REREFERENCING ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,CONFIG_STRUCT.fullRankAvRefBool);
%## saving EMG output
% tmp = pop_select(EEG,'channel',sort([eeg_chans,emg_chans]));
% tmp = pop_saveset(tmp,'filepath',mim_config.CleanEEG_output,...
%                     'filename',sprintf('%s_cleanEEG_%s_wEMG_wNoise',subj_char,preprocess_pipeline),...
%                     'savemode','twofiles');

%## save IMU, LS, & EEG output
% tmp = pop_select(EEG,'channel',sort([eeg_chans,biom_chans]));
% tmp = pop_saveset(tmp,'filepath',mim_config.CleanEEG_output,...
%                     'filename',sprintf('%s_clean_eeg_ls_imu.set',subj_char),...
%                     'savemode','twofiles');

%## SAVE JUST EEG OUTPUT
%(11/21/2024) Temporarily disabling this
EEG = pop_select(EEG,'channel',sort(eeg_chans));%2022-5-13 not use EMG
EEG = pop_saveset(EEG,'filepath',CONFIG_STRUCT.CleanEEG_output,...
                    'filename',sprintf('%s_cleanEEG_%s.set',subj_char,preprocess_pipeline),...
                    'savemode','twofiles');

%% HPG IMPLEMENTATION ONLY
% EEG = pop_select(EEG,'channel',sort([eeg_chans]));%2022-5-13 not use EMG
% disp([num2str(length(eeg_chans)),' remained']);
float_fPath = [EEG.filepath filesep sprintf('%s_cleanEEG_%s.fdt',subj_char,preprocess_pipeline)];
% [EEG,amica_cmd] = mim_prep_hpg_amica(EEG, ...
%     float_fPath, ...
%     mim_config.amicaEEG_output, ...
%     mim_config.emailStr, ...
%     avg_ref_pca_reduction);
[EEG,amica_cmd] = mim_prep_hpg_amica(EEG, ...
    float_fPath, ...
    CONFIG_STRUCT.amicaEEG_output, ...
    CONFIG_STRUCT.emailStr);
%--
EEG.etc.ica_pca_calc = sum(eig(cov(double(EEG.data'))) > 1e-7);
EEG = pop_saveset(EEG,'filepath',CONFIG_STRUCT.CleanEEG_output,...
                    'filename',sprintf('%s_cleanEEG_%s.set',subj_char,preprocess_pipeline),...
                    'savemode','twofiles');

%(11/21/2024) JS, Temporarily disabling this for imu/ls analysis
% EEG = pop_saveset(EEG,'filepath',mim_config.CleanEEG_output,...
%                     'filename',sprintf('%s_cleanEEG_%s.set',subj_char,preprocess_pipeline),...
%                     'savemode','twofiles');
toc
end

