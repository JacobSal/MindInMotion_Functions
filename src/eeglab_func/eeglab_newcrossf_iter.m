function [coh,mcoh,ttimes,lfreqs,cohboot,cohangles] = ...
    eeglab_newcrossf_iter(EEG1,EEG2,varargin)
%ADD_ANATOMICAL_LABELS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer:  Jacob Salminen
% Code Date: 04/30/2025, MATLAB 2023b
% Copyright (C) Jacob Salminen

%## TIME
tic
%## DEFINE DEFAULTS
DO_PRINTS = false;
%-- find fieldtrip on path
DEF_NEWCROSSF_STRUCT = struct(...
    'timewindow',[], ...
    'ntimesout',floor(EEG1.srate/pi()), ...
    'pointrange',1:size(EEG1.icaact,2), ...
    'baseline',nan(), ...
    'timelimits',[EEG1(1).xmin EEG1(1).xmax]*1000, ...
    'freqs',[3,250], ...
    'freqscale','log', ...
    'itctype','coher', ... % {'coher','phasecoher','phasecoher2'} 
    'do_plot_ersp',false, ...
    'do_plot_itc',false, ...
    'do_plot_phase',false, ...
    'do_single_dat',true, ...
    'do_timewarp',true, ...
    'cycles',[3,0.8], ...
    'padratio',4, ...
    'nfreqs',200, ...
    'maxfreq',250, ...
    'subitc','on', ...
    'alpha',nan(), ...
    'srate',EEG1.srate ...
    );
% NEWCROSSF_PKEEP = {'baseline','freqs','freqscale', ...
%     'padratio','nfreqs','alpha'};
NEWCROSSF_PKEEP = {'baseline','maxfreq','freqscale', ...
    'padratio','alpha','subitc'};

%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG1',@isstruct);
addRequired(p,'EEG2',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'NEWCROSSF_STRUCT',DEF_NEWCROSSF_STRUCT,@(x) validate_struct(x,DEF_NEWCROSSF_STRUCT,DO_PRINTS));
parse(p,EEG1,EEG2,varargin{:});
%## SET DEFAULTS
%-- extract vars
NEWCROSSF_STRUCT = p.Results.NEWCROSSF_STRUCT;
%-- set defaults
NEWCROSSF_STRUCT = set_defaults_struct(NEWCROSSF_STRUCT,DEF_NEWCROSSF_STRUCT,DO_PRINTS);

%% ===================================================================== %%
% fprintf('%s) Computing time/frequency decomposition...\n',EEG.subject);

%## COMPATABILITY EDITS FOR NEWTIMEF
onoffc = {'off','on'};
%--
newtimef_params = struct2args(NEWCROSSF_STRUCT);
fns = fieldnames(NEWCROSSF_STRUCT);
inds = find(cellfun(@(x) any(strcmp(x,NEWCROSSF_PKEEP)),fns));
keep_inds = zeros(length(newtimef_params),1);
for ff = 1:length(inds)
    ind = find(strcmp(fns(inds(ff)),newtimef_params));
    keep_inds(ind) = 1;
    keep_inds(ind+1) = 1;
end
newtimef_params = newtimef_params(logical(keep_inds));
% newtimef_params = {newtimef_params{:}; %#ok<CCAT>

%## CALC. TF & ITC
%--
% fprintf('processing component %i\n',k);

%-- make data compatible w/ newtimef and choose component
% tmpd = EEG.icaact(k,:,:);
tmpd1 = EEG1.icaact(1,NEWCROSSF_STRUCT.pointrange,:);
tmpd2 = EEG2.icaact(1,NEWCROSSF_STRUCT.pointrange,:);
tfdat1 = reshape(tmpd1,[1,length(NEWCROSSF_STRUCT.pointrange)*size(tmpd1,3)]);
tfdat2 = reshape(tmpd2,[1,length(NEWCROSSF_STRUCT.pointrange)*size(tmpd2,3)]);
tmpd1 = double.empty; %#ok<NASGU> % clear mem
tmpd2 = double.empty; %#ok<NASGU> % clear mem
EEG.icaact = double.empty; % clear mem

%## RUN NEWCROSSF
[coh,mcoh,ttimes,lfreqs,cohboot,cohangles] = ...
    newcrossf(tfdat1,tfdat2,...
    length(NEWCROSSF_STRUCT.pointrange), ...
    NEWCROSSF_STRUCT.timelimits, ...
    NEWCROSSF_STRUCT.srate, ...
    NEWCROSSF_STRUCT.cycles, ...
    newtimef_params{:});

if NEWCROSSF_STRUCT.do_single_dat
    coh = single(coh);
    mcoh = single(mcoh);
    cohboot = single(cohboot);
end

end

%% SUBFUNCTIONS ======================================================== %%
