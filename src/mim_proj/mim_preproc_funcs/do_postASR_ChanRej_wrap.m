function [EEG,EEG_temp_clean_timerej] = do_postASR_ChanRej_wrap(EEG,subjStr,cleaningMethod,do_postASR_ChanRej,MiM_config,fullRankAvRefBool)
        EEGpreprocess_fileName = [subjStr,'_cleanEEG_',do_postASR_ChanRej{1},'.set'];
        EEGpreprocess_filePath = fullfile(MiM_config.CleanEEG_output,do_postASR_ChanRej{1},subjStr);
        EEG_preprocessed = pop_loadset('filepath',EEGpreprocess_filePath,'filename',EEGpreprocess_fileName);
        [EEG,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,MiM_config);
        EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
        EEG_temp_clean_timerej = rerefC2CN2NExt2Ext_func(EEG_temp_clean_timerej,fullRankAvRefBool);
        preprocess_pipeline = [cleaningMethod];    
        mkdir(fullfile(MiM_config.CleanEEG_output,preprocess_pipeline,subjStr))
        fid = fopen(fullfile(MiM_config.CleanEEG_output,preprocess_pipeline,subjStr,'info.txt'),'w');
        fprintf(fid,['\n %.2f percent of frames were rejected\n'], p_frames_rej);
        fprintf(fid,['\n %.2f channels were rejected\n'], p_chan_rej);
        fclose(fid);
        zthreshold = EEG_temp_clean_timerej.etc.clean_artifacts.window_zthreshold;
        save(fullfile(MiM_config.CleanEEG_output,preprocess_pipeline,subjStr,'rejection_info.mat'),'p_frames_rej','p_chan_rej','zthreshold');

end