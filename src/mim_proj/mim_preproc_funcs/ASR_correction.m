function [EEG] = ASR_correction(EEG,params)
%MYASR Summary of this function goes here
%   Detailed explanation goes here
%   IN:
%       PATHS, 
%
%       ALLEEG,  
%
%       EEG, 
%           
%       CURRENTSET, 
%           
%       SUBJSTRUCT, 
%           
%       idx, 
%           
%       trialName, 
%            
%   OUT: 
%       ALLEEG,  
%
%       EEG, 
%           
%       CURRENTSET, 
%           
%   IMPORTANT: 

% Editted this file so that it works with merged file - Chang Liu
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {20,...
            false};
% p = inputParser;
%## REQUIRED
% addRequired(p,'EEG',@isstruct)
% addRequired(p,'fileList')
% %## OPTIONAL
% addOptional(p,'burstCritASR',Defaults{1})
% addOptional(p,'useExternalCalibASR',Defaults{2})
% parse(p, EEG, fileList, varargin{:});
%## SET DEFAULTS
burstCritASR = params.burstCritASR;
useExternalCalibASR = params.useExternalCalibASR;
FULL_RANK_AV_REF = false; %default false for regular avg ref, true for full rank version
%% --------------------------------------------------------------------- %%
myASRtic = tic;
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %redefine channels
EEG_only = pop_select(EEG,'channel',EEG_chans);
% restIndex = find(strcmpi('Rest.set',{fileList.name}));
StartEvents= find(strcmpi('TrialStart',{EEG.event.type}) & strcmpi('rest',{EEG.event.cond}));
EndEvents = find(strcmpi('TrialEnd',{EEG.event.type})& strcmpi('rest',{EEG.event.cond}));
startTime = EEG.times( round(EEG.event(StartEvents).latency) )/1000;
endTime = EEG.times( round(EEG.event(EndEvents).latency) )/1000;
disp(['Trial duration = ',num2str(endTime-startTime),' seconds']);
restEEG = pop_select( EEG, 'time',[startTime endTime ] );

%## if providing an external source of clean data. added check for rest existence
if useExternalCalibASR 
    cleanEEGforCalib_noChRej = restEEG; %ROEHL! need to make this smarter in the future. This is not dummy proof

%     cleanEEGforCalib_noChRej = pop_loadset('filename','NMM10_Clean_1.set','filepath',EEGfilePath); %load example of what clean EEG looks like
    cleanEEGforCalib = pop_select(cleanEEGforCalib_noChRej,'channel',{EEG_only.chanlocs.labels}); %ch rej the raw clean dataset to match the dataset you want to clean
    cleanEEGforCalib = rerefC2CN2NExt2Ext_func(cleanEEGforCalib,FULL_RANK_AV_REF); %re-reference the ch rejected clean dataset

    %NOTE: could do automatic time window rejection here (rather than
    %depend on user manual crop with excel sheet (Trial_Cropping.xls)
    %to help get the best clean resting data possible to use for our
    %calibration dataset. A subfunction of ASR does window rejection.
    %We can use that to clean up the resting dataset before becoming
    %calib data.

    %note you can use pop_clean_rawdata or clean_artifacts but pop_clean_rawdata has an error sometimes when you use external calibration data (error at line 159: com = sprintf('EEG = clean_artifacts(EEG, %s);', vararg2str(options));)
    EEG_ASR = clean_artifacts(EEG_only,...
        'BurstCriterionRefMaxBadChns',cleanEEGforCalib,...
        'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off',...
        'Highpass','off','BurstCriterion',burstCritASR,'WindowCriterion','off',...
        'BurstRejection','off','Distance','Euclidian');
else
    %let ASR find calibration on its own (taken from same file as the one you are trying to clean)
    %run asr at desired burst setting and keep track of time
    EEG_ASR = clean_artifacts(EEG, ...
        'FlatlineCriterion','off', 'ChannelCriterion','off','LineNoiseCriterion','off',...
        'Highpass','off','BurstCriterion',burstCritASR,'WindowCriterion','off',...
        'BurstRejection','off','Distance','Euclidian');

end

% vis_artifacts(EEG_ASR,EEG_only);
EEG.data(EEG_chans,:) = EEG_ASR.data; %copy over cleaned data (note: can't handle ch rej with ASR right now)
clear EEG_ASR;

%## Clean up EEG for AMICA
EEG = rmfield(EEG,'icaweights');
EEG.icaweights = [];
EEG = rmfield(EEG,'icawinv');
EEG.icawinv = [];
EEG = rmfield(EEG,'icaact');
EEG.icaact = [];
EEG = rmfield(EEG,'icasphere');
EEG.icasphere = [];
EEG = rmfield(EEG,'icachansind');
EEG.icachansind = [];%also remove icachanind?

executionTime = toc(myASRtic);
disp(['ASR cleaning took ',num2str(executionTime),' seconds']);
end

