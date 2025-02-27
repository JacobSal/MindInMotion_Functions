function [EEG_out,Total_rej] = autoRejCh_func_CL(EEG,thres,EEG_chans,Noise_chans,EMG_chans)
% Editted by Chang - add ability to specify threshold
%% define channels
% All_chans = 1:EEG.nbchan;
% EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
% Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
% EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
% BioM_chans = find(strcmpi('BioM',{EEG.chanlocs.type}));
% Acc_chans = find(strcmpi('Acc',{EEG.chanlocs.type}));
% EOG_chans = find(strcmpi('EOG',{EEG.chanlocs.type}));
% NEOG_chans = find(strcmpi('N-EOG',{EEG.chanlocs.type}));
% Sync_chans = find(strcmpi('Sync',{EEG.chanlocs.type}));
% Other_chans = setdiff(All_chans, unique([EEG_chans, Noise_chans, EMG_chans, BioM_chans, Acc_chans, EOG_chans, NEOG_chans, Sync_chans]));


%% Identify bad channels within each channel type of interest
% badEEGch = EEG_chans( find(EEGchRange > 500 | EEGchRange < 5) )
% badEMGch = EMG_chans( find(EMGchRange > 3000 | EMGchRange < 5) )
% badNoiseCh = Noise_chans( find(NoiseChRange > 1000 | NoiseChRange < 5) )


%new standard deviation method
EEG_SD = std(EEG.data(EEG_chans,:)');
EMG_SD = std(EEG.data(EMG_chans,:)');
Noise_SD = std(EEG.data(Noise_chans,:)');

badEEGch = EEG_chans( find(EEG_SD > thres*median(EEG_SD)) )
badEMGch = EMG_chans( find(EMG_SD > thres*median(EMG_SD)) )
badNoiseCh = Noise_chans( find(Noise_SD > thres*median(Noise_SD)) )

%% Display
fprintf('We are rejecting EEG channels %s \n',num2str(badEEGch));
fprintf('We are rejecting EMG channels %s \n',num2str(badEMGch));
fprintf('We are rejecting Noise channels %s \n',num2str(badNoiseCh));

%% reject chans
% if ~isfield(EEG.etc,'remove_chan')
%    EEG.etc.remove_chan = [];
%    EEG.etc.remove_chan = EEG.chanlocs([badEEGch badEMGch badNoiseCh]);
% else
%    EEG.etc.remove_chan = vertcat(EEG.etc.remove_chan,EEG.chanlocs([badEEGch badEMGch badNoiseCh]));
% end

%(11/21/2024) JS, Keeping IMU & LS data
% EEG_out = pop_select( EEG,'nochannel',sort( unique([badEEGch, badEMGch, badNoiseCh]) ) );
EEG_out = pop_select( EEG,'nochannel',sort( unique([badEEGch, badEMGch, badNoiseCh]) ) );

%% Output the number of channels got rejected
Total_rej = [length(badEEGch) length(badEMGch) length(badNoiseCh)];
end