function [ersp_dat,itc_dat,mbase,ttimes,lfreqs,ersp_boot,itc_boot,trial_ersp_dat] = ...
    eeglab_newtimef_iter(EEG,k,varargin)
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
DEF_NEWTIMEF_STRUCT = struct(...
    'timewindow',[], ...
    'timewarp',EEG.timewarp.latencies, ...
    'timewarpms',EEG.timewarp.warpto, ...
    'ntimesout',floor(EEG.srate/pi()), ...
    'pointrange',EEG.timewarp.epochs, ...
    'baseline',nan(), ...
    'timelimits',[EEG(1).xmin EEG(1).xmax]*1000, ...
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
    'alpha',nan(), ...
    'srate',EEG.srate, ...
    'newtimef_params',{{}} ...
    );
NEWTIMEF_TWPKEEP = {'timewarp','timewarpms', ...
    'baseline','freqs','freqscale', ...
    'padratio','nfreqs','alpha','ntimesout','itctype'};
NEWTIMEF_RPKEEP = {'baseline','freqs','freqscale', ...
    'padratio','nfreqs','alpha','ntimesout','itctype'};

%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'k',@isnumeric)
%## OPTIONAL
%## PARAMETER
addParameter(p,'NEWTIMEF_STRUCT',DEF_NEWTIMEF_STRUCT,@(x) validate_struct(x,DEF_NEWTIMEF_STRUCT,DO_PRINTS));
parse(p,EEG,k,varargin{:});
%## SET DEFAULTS
%-- extract vars
NEWTIMEF_STRUCT = p.Results.NEWTIMEF_STRUCT;
%-- set defaults
NEWTIMEF_STRUCT = set_defaults_struct(NEWTIMEF_STRUCT,DEF_NEWTIMEF_STRUCT,DO_PRINTS);

%% ===================================================================== %%
% fprintf('%s) Computing time/frequency decomposition...\n',EEG.subject);

%## COMPATABILITY EDITS FOR NEWTIMEF
onoffc = {'off','on'};
%--
newtimef_params = struct2args(NEWTIMEF_STRUCT);
fns = fieldnames(NEWTIMEF_STRUCT);
if NEWTIMEF_STRUCT.do_timewarp
    inds = find(cellfun(@(x) any(strcmp(x,NEWTIMEF_TWPKEEP)),fns));
    keep_inds = zeros(length(newtimef_params),1);
    for ff = 1:length(inds)
        ind = find(strcmp(fns(inds(ff)),newtimef_params));
        keep_inds(ind) = 1;
        keep_inds(ind+1) = 1;
    end
    newtimef_params = newtimef_params(logical(keep_inds));
else
    inds = find(cellfun(@(x) any(strcmp(x,NEWTIMEF_RPKEEP)),fns));
    keep_inds = zeros(length(newtimef_params),1);
    for ff = 1:length(inds)
        ind = find(strcmp(fns(inds(ff)),newtimef_params));
        keep_inds(ind) = 1;
        keep_inds(ind+1) = 1;
    end
    newtimef_params = newtimef_params(logical(keep_inds));
end 
newtimef_params = {newtimef_params{:}, ...
    'plotersp', onoffc{double(NEWTIMEF_STRUCT.do_plot_ersp)+1}, ...
    'plotitc', onoffc{double(NEWTIMEF_STRUCT.do_plot_itc)+1}, ...
    'plotphase', onoffc{double(NEWTIMEF_STRUCT.do_plot_phase)+1}, ...
    'verbose', 'off', ...
    NEWTIMEF_STRUCT.newtimef_params{:}}; %#ok<CCAT>
%## CALC. TF & ITC
%--
% fprintf('processing component %i\n',k);

%-- make data compatible w/ newtimef and choose component
tmpd = EEG.icaact(:,NEWTIMEF_STRUCT.pointrange,:);
timefdata  = reshape(tmpd,[1,length(NEWTIMEF_STRUCT.pointrange)*size(tmpd,3)]);
tmpd = double.empty; %#ok<NASGU> % clear mem
EEG.icaact = double.empty; % clear mem

%## RUN NEWTIMEF
if ~isempty(timefdata)
    [ersp_dat,itc_dat,mbase,ttimes,lfreqs,ersp_boot,itc_boot,trial_ersp_dat,~] ...
          = newtimef(timefdata, ...
          length(NEWTIMEF_STRUCT.pointrange), ...
          NEWTIMEF_STRUCT.timelimits, ...
          NEWTIMEF_STRUCT.srate, ...
          NEWTIMEF_STRUCT.cycles, ...
          newtimef_params{:});
end

if NEWTIMEF_STRUCT.do_single_dat
    trial_ersp_dat = single(trial_ersp_dat);
    ersp_dat = single(ersp_dat);
    itc_dat = single(itc_dat);
    ersp_boot = single(ersp_boot);
    itc_boot = single(itc_boot);
end

end

%% SUBFUNCTIONS ======================================================== %%
