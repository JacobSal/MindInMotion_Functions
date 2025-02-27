function [ersp_tf,freqs,trialinfo,ersp_nolog_tf] = spca_get_ersp_dat(tf_mat_struct,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## Define Parser
condition_char = [];
p = inputParser;
%## REQUIRED
addRequired(p,'f_mat_struct',@isstruct);
%## OPTIONAL
addOptional(p,'condition_char',condition_char,@ischar);
%## PARAMETER
parse(p,tf_mat_struct,varargin{:});
%## SET DEFAULTS
condition_char = p.Results.condition_char;
%% (ALTERNATIVE) USE TIMEWARPPING EEGLAB =================================== %%
% functions: newtimef.m, timefreq.m, timewarp.m
% (12/9/2023) JS, the timewarp function is the real meat that is needed to
% warp gait events to a grand average. Could use the morlet_transform_fast
% after that.
%## GET ERSP DATA FROM ICASPEC
%- reshape data [trials x pnts x chans x freq]
fn = fieldnames(tf_mat_struct);
inds = find(contains(fieldnames(tf_mat_struct),'comp'));
test = tf_mat_struct.(fn{inds(1)});
%- store in matrix
ersp_nolog_tf = zeros(size(test,1),size(test,2),size(test,3),length(inds),'double');
for i = 1:length(inds)
    ersp_nolog_tf(:,:,:,i) = tf_mat_struct.(fn{inds(i)}); % freq x epoch x chan
end
%- extract specific condition
if ~isempty(condition_char)
    inds = strcmp({tf_mat_struct.trialinfo.cond},condition_char);
    ersp_nolog_tf = squeeze(ersp_nolog_tf(:,:,inds,:));
end
%- extract other useful information
freqs = tf_mat_struct.freqs;
trialinfo = tf_mat_struct.trialinfo;
ersp_tf = 20*log10(ersp_nolog_tf);
end

