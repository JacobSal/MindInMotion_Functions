
function [clean_data, clean_data_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,params)
%MAIN_FUNC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%
% NOTES:
%       tmp_eeg_clean: is with only channel rejection
%       tmp_eeg_clean_timerej:  with both channel rej and window rejection
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/23/2023, MATLAB 2019a
% Copyright (C) Noelle Jacobsen, Designer & Creator
% Copyright (C) Chang Liu, Designer & Creator (20210820)
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu, Code Dweeb
% 
% Modified - 20230223 - function wrapper (Jacob Salminen)
plot_chan_rej = 0;
%% Reject bad EEG channels and time rejection
%(11/21/2024) JS, originally commented out
%[clean_data, clean_artifact_data] = channelrejection(original_data, Options...)
% Detect bad channels in EEG
%## DEFINE CURRENT CHANNELS
eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %redefine channels
noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
emg_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
biom_chans = find(strcmpi('BioM',{EEG.chanlocs.type}));

eeg_only_set = pop_select(EEG,'channel',sort(eeg_chans));
emg_noise_biom_set = pop_select(EEG,'channel',sort([emg_chans,noise_chans,biom_chans]));
% biom_set = pop_select(EEG,'channel',sort([biom_chans]));

%## CHANNEL REJECTION & TIME REJECTION
%- use the channel rejection function from Noelle
[tmp_eeg_clean,tmp_eeg_clean_timerej] = channelrejection(eeg_only_set,...
    'burst_crit','off',...
    'wind_crit',params.wind_crit,...
    'kurt_crit',params.kurt_crit,...
    'chan_crit1',params.chan_crit1,...
    'flatline_crit1',params.flatline_crit,...
    'highpass',[0.25 0.75],'flatline_crit',5,...
    'line_crit','off','channel_crit_maxbad_time',0,...
    'chan_range',find(strcmpi('EEG',{eeg_only_set.chanlocs.type})),...
    'chan_detect_num_iter',1,...%params.chan_detect_num_iter
    'window_crit_tolerances',params.window_crit_tolerances,...% Used to be [-inf,10] 20220307
    'kurt_prob_std_rm',params.kurt_prob_std_rm,...
    'min_clean_fraction',0.25);
disp('Done with clean artifacts');

%- Compute percent of frames rejected
frames_before = length(eeg_only_set.times);
frames_after = length(tmp_eeg_clean_timerej.times);
p_frames_kept = (frames_after/frames_before)*100;
p_frames_rej = 100-p_frames_kept;
EEG.etc.percent_frames_rej= p_frames_rej;
fprintf('\n %.2f percent of frames were rejected\n', p_frames_rej);

%- Store params in etc file 
fprintf('Saving the cleaned set')
tmp_eeg_clean.etc.notes = params; % SAVE the pre-processing parameters
tmp_eeg_clean.etc.RejTimeWindow = 0;
tmp_eeg_clean_timerej.etc.notes = params;
tmp_eeg_clean_timerej.etc.RejTimeWindow = 1;

%     [ALLEEG EEG] = eeg_store(ALLEEG,EEG_temp_clean,0);
%     [ALLEEG EEG] = eeg_store(ALLEEG,EEG_temp_clean_timerej,0);
numChan_before = length(eeg_only_set.chanlocs);
numChan_after = length(tmp_eeg_clean_timerej.chanlocs);
p_chan_rej = numChan_after - numChan_before;
fprintf('\n %.2f channels were rejected\n', p_chan_rej);
%- redefine chans
eeg_chans = find(strcmpi('EEG',{tmp_eeg_clean_timerej.chanlocs.type})); %redefine channels

%% Combine the cleaned EEG and with the noise/EMG/BioM data
%- Add data into the EEG only set
try
    %- get clean EEG data & combine it with the EMG, BioM, and Noise data
    tmp_eeg_clean_timerej.data = [tmp_eeg_clean_timerej.data;
        emg_noise_biom_set.data(:,find(tmp_eeg_clean_timerej.etc.clean_artifacts.clean_sample_mask == 1))];
    
catch
    %- catch a bug where the rejection is one sample off
    disp('omitting the last sample')
    clean_indx = find(tmp_eeg_clean_timerej.etc.clean_artifacts.clean_sample_mask == 1);
    tmp_eeg_clean_timerej.data = [tmp_eeg_clean_timerej.data;
        emg_noise_biom_set.data(:,clean_indx(1:end-1))];
end
%- Add EMG, Noise, & BIOM data into non-timerej set.
tmp_eeg_clean.data = [tmp_eeg_clean.data; emg_noise_biom_set.data];    
n_eeg_chans = length(eeg_chans);
for chan = 1:length(sort([emg_chans,noise_chans,biom_chans]))
    %- fix not consistent chan info
    tmp_eeg_clean_timerej.chanlocs(n_eeg_chans+chan).urchan = [];
    tmp_eeg_clean_timerej.chanlocs(n_eeg_chans+chan).ref = [];
    tmp_eeg_clean_timerej.chanlocs(n_eeg_chans+chan).sph_theta_besa = [];
    tmp_eeg_clean_timerej.chanlocs(n_eeg_chans+chan).sph_phi_besa = [];

    tmp_eeg_clean.chanlocs(n_eeg_chans+chan).urchan = [];
    tmp_eeg_clean.chanlocs(n_eeg_chans+chan).ref = [];
    tmp_eeg_clean.chanlocs(n_eeg_chans+chan).sph_theta_besa = [];
    tmp_eeg_clean.chanlocs(n_eeg_chans+chan).sph_phi_besa = [];

    if ~isfield(emg_noise_biom_set.chanlocs(chan),'sph_theta_besa')
        emg_noise_biom_set.chanlocs(chan).sph_theta_besa = [];
        emg_noise_biom_set.chanlocs(chan).sph_phi_besa = [];
    end
    tmp_eeg_clean_timerej.chanlocs(n_eeg_chans+chan) = emg_noise_biom_set.chanlocs(chan);
    tmp_eeg_clean.chanlocs(n_eeg_chans+chan) = emg_noise_biom_set.chanlocs(chan);
end
% update other information as well
tmp_eeg_clean_timerej.nbchan = size(tmp_eeg_clean_timerej.data,1);
tmp_eeg_clean_timerej.comments = pop_comments(EEG.comments,'','Cleaned EEG + Noise + EMG + BioM');

%% Visualization after channel rejection to check if the channels removed make sense
% visualize the cleaning, in order to plot the channels we removed, we must
% use the clean_channel_mask function
% clean_chann
origChannelInfo = cellstr(strvcat(EEG.urchanlocs.labels));
update_Channels = cellstr(strvcat(tmp_eeg_clean_timerej.chanlocs.labels));
% idx_bad_chans = cellfun(@(s) find(strcmpi(s, origChannelInfo)), ...
%     update_Channels,'UniformOutput', true);
clean_channel_mask = zeros(length(origChannelInfo),1);
if length(update_Channels) ~= length(origChannelInfo)
    idx_chans_to_remove = cellfun(@(s) find(strcmpi(s, origChannelInfo)), update_Channels,'UniformOutput', true);
else
    idx_chans_to_remove = [];
end
%- (11/21/2024) JS, this code doesn't seem to do anything
clean_channel_mask(idx_chans_to_remove) = 1;
tmp_eeg_clean.etc.channel_mask = logical(clean_channel_mask);

tmp_eeg_clean_timerej.etc.channel_mask = logical(clean_channel_mask);
tmp_eeg_clean_timerej.etc.clean_sample_mask = logical(tmp_eeg_clean_timerej.etc.clean_artifacts.clean_sample_mask);

EEG.etc.clean_channel_mask = logical(clean_channel_mask);
EEG.etc.clean_sample_mask = logical(tmp_eeg_clean_timerej.etc.clean_artifacts.clean_sample_mask);

%## (11/21/2024) JS, WARNING this is not working
%{
if plot_chan_rej
    % NEED to fix
    EEG_chans_plot = find(strcmpi('EEG',{EEG.urchanlocs.type})); %redefine channels (numbering changed since we took subset of channels)

    EEG = rerefC2CN2NExt2Ext_func(EEG,params.fullRankAvRefBool);
    EEG_orig_reref = rerefC2CN2NExt2Ext_func(EEG_orig,params.fullRankAvRefBool);
    disp('plotting now')
    % if getting stuck, remove eeglab path from R drive. % Wont work if
    % there is spare amp
    vis_artifacts(EEG,EEG_orig_reref,'show_events',0,...
        'ChannelSubset',find(strcmpi('EEG',{EEG_orig_reref.chanlocs.type})),...
        'show_removed_portions',1,'EqualizeChannelScaling',1);
    
    %test
    % TO DO ADD PLOTSSSS, and Iterative rejection to avoid rejecting 'okay'
    % components (BeMoBi pipeline)
    bemobil_plot_EEG_visartifact(EEG,EEG_orig_reref)
    saveas(gcf,strcat(fullfile(params.MIMDataFolder_CL,subjStr,'EEG','Figures'),'\pre-post processing EEG_',t,'.jpg'));

    % TO DO: add spectra plots
end
%}

%% SAVE after channel rejection before proceeding
% TO DO: 
%     EEG_chans = find(strcmpi('EEG',{EEG_temp_clean_timerej.chanlocs.type})); %redefine channels (numbering changed since we took subset of channels)
%     EMG_chans = find(strcmpi('EMG',{EEG_temp_clean_timerej.chanlocs.type}));
%     Noise_chans = find(strcmpi('Noise',{EEG_temp_clean_timerej.chanlocs.type}));
% 
%     %update channel labels again just in case?
%     EEG_chans = find(strcmpi('EEG',{EEG_temp_clean.chanlocs.type})); %redefine channels (numbering changed since we took subset of channels)
%     EMG_chans = find(strcmpi('EMG',{EEG_temp_clean.chanlocs.type}));
%     Noise_chans = find(strcmpi('Noise',{EEG_temp_clean.chanlocs.type}));

tmp_eeg_clean_timerej.comment = 'clean_artifact_window_rejection';
tmp_eeg_clean_timerej = rerefC2CN2NExt2Ext_func(tmp_eeg_clean_timerej,params.fullRankAvRefBool);

tmp_eeg_clean.comment = 'clean_artifact_channel_rejection';
tmp_eeg_clean = rerefC2CN2NExt2Ext_func(tmp_eeg_clean,params.fullRankAvRefBool);

clean_data = tmp_eeg_clean;
clean_data_timerej = tmp_eeg_clean_timerej;

end
