function [EEG] = mim_load_one_subj(fname,fpath,subj_char,cond_char,group_char,sess_char,varargin)
%MIM_CREATE_ALLEEG Summary of this function goes here
%   Detailed explanation goes here
%     MATLAB CODE DOCUMENTATION

%     Code Designer: Jacob Salminen
%     Date: 11/25/2022
% 
%     Functions:
%     1. fem_eeglab_dipfit: Configures dipole fitting parameters for Finite Element Model (FEM) dipfit analysis.
%     2. custom_update_chanlocs: Updates EEG channel locations using custom electrode locations.
%     3. bem_eeglab_dipfit: Fits `dipole models to EEG data using Boundary Element Model (BEM) dipfit analysis.
% 
%     See individual function descriptions for more details on their functionality.
%- logo
% cat_logo();
%- time
tt = tic;
%## PATHS
if ~ispc
    pp = strsplit(path,':');
else
    pp = strsplit(path,';');
end
%-
tmp = regexp(pp,'eeglab');
tmp = pp(~cellfun(@isempty,tmp));
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_eeglab = b2{1}(1:b1{1});
    fprintf('EEGLAB path: %s\n',path_eeglab);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('EEGLAB path not found.\n');
    end
end
%## TIME
tic
%## DEFINE DEFAULTS
DEF_ERSP_STRUCT = struct('components',[], ...
    'channels',{}, ...
    'powbase',[], ...
    'trialindices',[], ...
    'savetrials',{'off'}, ... % 'on' | 'off'
    'plot',{'off'}, ... % 'on' | 'off'  %NOTE: not documented for debugging purpose
    'recompute',{'off'}, ... % 'on' | 'off'
    'getparams',{'off'}, ... % 'on' | 'off'
    'savefile',{'on'}, ... % 'on' | 'off'
    'parallel',{'off'}, ... % 'on' | 'off'
    'timewindow',[], ... %NOTE: ignored, deprecated
    'fileout', {''}, ...
    'timelimits',[EEG(1).xmin EEG(1).xmax]*1000, ...
    'cycles',[3 .5], ...
    'padratio',1, ...
    'trialinfo',struct([]), ...
    'freqs',[0 EEG(1).srate/2], ...
    'rmcomps',{cell(1,length(EEG))}, ...
    'interp',struct([]), ...
    'freqscale',{'log'}, ...
    'alpha',NaN, ...
    'baseline',0, ...
    'type',{'ersp'}); % 'ersp','itc','both','ersp&itc'

%##
DEF_NEWTIMEF_PARAMS = struct('boottype','shuffle', ... % {'shuffle','rand','randall'}  
    'condboot','abs', ... % {'abs','angle','complex'}
    'title',{{'ERSP'}}, ...
    'title2',{'ERSP'}, ...
    'winsize',DEFAULT_WINSIZE, ... % [0 Inf]
    'pad',DEFAULT_PAD, ...
    'timesout',DEFAULT_NWIN, ...
    'padratio',DEFAULT_OVERSMP, ... % [0 Inf]
    'topovec',[], ...
    'elocs',DEFAULT_ELOC, ... % {'string','struct'}
    'alpha',DEFAULT_ALPHA, ... % [0 0.5]
    'marktimes',DEFAULT_MARKTIME, ...
    'powbase',NaN, ...
    'pboot',NaN, ...
    'rboot',NaN, ...
    'plotersp','on', ... % 'on' | 'off'
    'plotamp','on', ... % 'on' | 'off'
    'plotitc','on', ... % 'on' | 'off'
    'detrend','off', ... % 'on' | 'off'
    'rmerp','off', ... % 'on' | 'off'
    'basenorm','off', ... % 'on' | 'off'
    'commonbase','on', ... % 'on' | 'off'
    'baseline',0, ...
    'baseboot',1, ...
    'linewidth',2, ...
    'naccu',200, ... % [1 Inf]
    'mtaper',[], ...
    'maxfreq',DEFAULT_MAXFREQ, ... % [0 Inf]
    'freqs',[0 DEFAULT_MAXFREQ], ... % [0 Inf]
    'cycles',[], ...
    'nfreqs',[], ...
    'freqscale','linear', ...
    'vert',[],  ...
    'newfig','on', ... % 'on' | 'off'
    'type','phasecoher', ... % {'coher','phasecoher','phasecoher2'}
    'itctype','phasecoher', ... % {'coher','phasecoher','phasecoher2'}
    'outputformat','plot', ... % {'old','new','plot' }
    'phsamp','off', ...  % 'on' | 'off' % NOTE: phsamp not completed - Toby 9.28.2006
    'plotphaseonly','off', ... % 'on' | 'off'
    'plotphasesign','on', ... % 'on' | 'off'
    'plotphase','on', ... % 'on' | 'off' % NOTE: same as above for backward compatibility
    'pcontour','off', ... % 'on' | 'off'
    'precomputed',struct([]), ...
    'itcmax',[], ...
    'erspmax',[], ...
    'lowmem','off', ... % 'on' | 'off'
    'verbose','on', ... % 'on' | 'off'
    'plottype','image', ... % {'image','curve'}
    'mcorrect','none', ... % {'fdr','none'}
    'plotmean','on', ... % 'on' | 'off'
    'plotmode','', ... % for metaplottopo
    'highlightmode','background', ... % {'background','bottom'}
    'chaninfo',struct([]), ...
    'erspmarglim',double.empty, ...
    'itcavglim',double.empty, ...
    'erplim',double.empty, ...
    'speclim',double.empty, ...
    'ntimesout',[], ...
    'scale','log', ... % { 'log','abs'}
    'timewarp',double.empty, ...
    'timewarpms',double.empty, ...
    'timewarpfr',double.empty, ...
    'timewarpidx',double.empty, ...
    'timewarpidx',double.empty, ...
    'timeStretchMarks',double.empty, ...
    'timeStretchRefs',double.empty, ...
    'timeStretchPlot',double.empty, ...
    'trialbase','off', ...  % 'on' | 'off' | 'full'
    'caption','', ...
    'hzdir',HZDIR, ... % {'up','down','normal','reverse'}
    'ydir',YDIR, ... % {'up','down','normal','reverse'}
    'cycleinc','linear', ... % {'linear','log'} 
    'colormap',DEFAULT_COLORMAP,...
    );
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## PARAMETER
addParameter(p,'ERSP_STRUCT',DEF_ERSP_STRUCT,@(x) validate_struct(x,DEF_ERSP_STRUCT));
parse(p,fname,fpath,subj_char,cond_char,group_char,sess_char,varargin{:});
%## SET DEFAULTS
ERSP_STRUCT = p.Results.ERSP_STRUCT;
%-
ERSP_STRUCT = set_defaults_struct(ERSP_STRUCT,DEF_ERSP_STRUCT);
%% ===================================================================== %%
tmpparams = parameters;
if length(g.indices) > 1
    tmpparams{end+1} = 'verbose';
    tmpparams{end+1} = 'off';
end

% Run timef() to get ERSP
% ------------------------
timefdata  = reshape(X(k,pointrange,:), 1, length(pointrange)*size(X,3));
mytimes = [];
mylogfreqs = [];
alltfX   = [];
if ~isempty(timefdata)
    [logersp,logitc,logbase,mytimes,mylogfreqs,logeboot,logiboot,alltfX] ...
          = newtimef( timefdata, length(pointrange), g.timelimits, EEG(1).srate, tmpparams{2:end});
end
%if strcmpi(g.plot, 'on'), return; end
if usesingle
    alltfX = single(alltfX);
end

allTrialsTmp{k}   = single( alltfX );
allTrialsTime{k}  = mytimes;
allTrialsFreqs{k} = mylogfreqs;
