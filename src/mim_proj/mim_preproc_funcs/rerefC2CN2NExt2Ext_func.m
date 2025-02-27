function EEG = rerefC2CN2NExt2Ext_func(EEG,full_rank_bool)
%re-reference cort to cort, noise to noise, and 

if full_rank_bool
    disp('Rereferencing cort to cort, noise to noise, and emg to emg (with ZERO loss of rank!)');
else
    disp('Rereferencing cort to cort, noise to noise, and emg to emg (with a loss of rank of 1 for each rereferencing)');
end

%% define channels
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type})); 
%note: may need to add more lines in future if other external sensors are
%used with a different channel type than 'EMG' (e.g. could use an array of
%EOG sensors and want to re-reference them to themselves)

%% re-ref cort to cort and noise to noise and emg to emg 
% EEG = pop_reref( EEG, [],'exclude',Noise_chans );
% EEG = pop_reref( EEG, [],'exclude',sort([EEG_chans, EMG_chans]) );
%(02/21/2025) JS, there is a bug in these functions that causes rank
%deficiency issues. See. 
% Kim, H., Luo, J., Chu, S., Cannard, C., Hoffmann, S., & Miyakoshi, M. 
% (2023). ICA’s bug: How ghost ICs emerge from effective rank deficiency 
% caused by EEG electrode interpolation and incorrect re-referencing. 
% Frontiers in Signal Processing, 3. https://doi.org/10.3389/frsip.2023.1064138

%## CALC REREFENCING
chansin = Noise_chans;
nchansin = length(chansin);
if full_rank_bool
    refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%new full rank method 
    %Note CL: creating a matrix so that EEG*refmatrix = new_EEG(after referencing)
else
    refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; % normal avg ref (rank losing method)
end
chansout = chansin;
EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
chansin = EEG_chans;
nchansin = length(chansin);
if full_rank_bool
    refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);% new full rank method 
else
    refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; % normal avg ref (rank losing method)
end
chansout = chansin;
EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
chansin = EMG_chans;
nchansin = length(chansin);
if full_rank_bool
    refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);% new full rank method 
else
    refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; % normal avg ref (rank losing method)
end
chansout = chansin;
EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
end