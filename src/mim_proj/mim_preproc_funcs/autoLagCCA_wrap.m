function EEG = autoLagCCA_wrap(EEG,params)
    getChannelTypes; 
    EEGandEMG = pop_select(EEG,'channel',sort([EEG_chans, EMG_chans]));%'channel' = retain, why retain EMG
    [EEG_CCA] = autoLagCCA(EEGandEMG,params.lagAmount_samples ,params.CCA_Rsq_thres);
    % autoLagCCA_figure(EEG_CCA)
    disp('Done CCA');
    EEG.data(sort([EEG_chans,EMG_chans]),:) = EEG_CCA.data; %copy over cleaned data  %pay attention that indices line up (EMG interspersed with EEG in terms of chan numbers)
    clear EEG_CCA; EEG.comments = 'CCA after iCC';
end