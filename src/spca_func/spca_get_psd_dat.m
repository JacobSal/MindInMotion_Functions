function [psd_f,freqs,trialinfo,psd_nolog_f] = spca_get_psd_dat(f_mat_struct,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## Define Parser
COND_STR = [];
p = inputParser;
%## REQUIRED
addRequired(p,'f_mat_struct',@isstruct);
addRequired(p,'base_f_mean',@isnumeric);
%## OPTIONAL
addParameter(p,'COND_STR',COND_STR,@ischar);
%## PARAMETER
parse(p,f_mat_struct,varargin{:});
%## SET DEFAULTS
COND_STR = p.Results.COND_STR;
%% (ALTERNATIVE) USE TIMEWARPPING EEGLAB =================================== %%
% functions: newtimef.m, timefreq.m, timewarp.m
% (12/9/2023) JS, the timewarp function is the real meat that is needed to
% warp gait events to a grand average. Could use the morlet_transform_fast
% after that.
%## GET PSD DATA FROM ICASPEC
%- reshape data [trials x pnts x chans x freq]
fn = fieldnames(f_mat_struct);
inds = find(contains(fieldnames(f_mat_struct),'comp'));
test = f_mat_struct.(fn{inds(1)});
%- store in matrix
psd_f = zeros(size(test,1),size(test,2),length(inds),'double');
for i = 1:length(inds)
    psd_f(:,:,i) = f_mat_struct.(fn{inds(i)}); % freq x epoch x chan
end
%- extract specific condition
if ~isempty(COND_STR)
    inds = strcmp({f_mat_struct.trialinfo.cond},COND_STR);
    psd_f = squeeze(psd_f(:,inds,:));
end
%- extract other useful information
freqs = f_mat_struct.freqs;
trialinfo = f_mat_struct.trialinfo;
psd_nolog_f = 10.^(eeg_psd/10);
end

