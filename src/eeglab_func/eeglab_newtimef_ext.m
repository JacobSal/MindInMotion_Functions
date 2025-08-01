function [timef_data_struct, timef_dat, times, logfreqs, params] = eeglab_newtimef_ext(STUDY,EEG,varargin)
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
%-- find fieldtrip on path
DEF_NEWTIMEF_STRUCT = struct(...
    'components',1:size(EEG.icaact,1), ...
    'timewindow',[], ...
    'timewarp',EEG.timewarp.latencies, ...
    'timewarpms',EEG.timewarp.warpto, ...
    'ntimesout',EEG.srate/pi(), ...
    'pointrange',EEG.timewarp.epochs, ...
    'baseline',nan(), ...
    'timelimits',[EEG(1).xmin EEG(1).xmax]*1000, ...
    'freqs',[3,250], ...
    'freqscale','log', ...
    'indices',1:size(EEG.icaact,1), ...
    'do_plot_ersp',false, ...
    'do_plot_itc',false, ...
    'do_plot_phase',false, ...
    'cycles',[3,0.8], ...
    'padratio',1, ...
    'nfreqs',200, ...
    'alpha',nan(), ...
    'srate',EEG.srate ...
    );
NEWTIMEF_PKEEP = {'timewarp','timewarpms', ...
    'baseline','timelimits','freqs','freqscale', ...
    'padratio','nfreqs','alpha','ntimesout'};
%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'NEWTIMEF_STRUCT',DEF_NEWTIMEF_STRUCT,@(x) validate_struct(x,DEF_NEWTIMEF_STRUCT));
parse(p,STUDY,EEG,varargin{:});
%## SET DEFAULTS
%-- extract vars
NEWTIMEF_STRUCT = p.Results.NEWTIMEF_STRUCT;
%-- set defaults
NEWTIMEF_STRUCT = set_defaults_struct(NEWTIMEF_STRUCT,DEF_NEWTIMEF_STRUCT);

%% ===================================================================== %%
fprintf('%s) Computing time/frequency decomposition...\n',EEG.subject);

%## COMPATABILITY EDITS FOR NEWTIMEF
onoffc = {'off','on'};
%--
newtimef_params = struct2args(NEWTIMEF_STRUCT);
fns = fieldnames(NEWTIMEF_STRUCT);
inds = find(cellfun(@(x) any(strcmp(x,NEWTIMEF_PKEEP)),fns));
keep_inds = zeros(length(newtimef_params),1);
for ff = 1:length(inds)
    ind = find(strcmp(fns(inds(ff)),newtimef_params));
    keep_inds(ind) = 1;
    keep_inds(ind+1) = 1;
end
newtimef_params = newtimef_params(logical(keep_inds));
newtimef_params = {newtimef_params{:}, ...
    'plotersp', onoffc{double(NEWTIMEF_STRUCT.do_plot_ersp)+1}, ...
    'plotitc', onoffc{double(NEWTIMEF_STRUCT.do_plot_itc)+1}, ...
    'plotphase', onoffc{double(NEWTIMEF_STRUCT.do_plot_phase)+1}, ...
    'verbose', 'off'}; %#ok<CCAT>

%## LOOP THROUGH COMPONENTS & CALC. TF
datin = EEG.icaact;
do_single_dat = true;
datout = cell(length(NEWTIMEF_STRUCT.components));
timesout = cell(length(NEWTIMEF_STRUCT.components));
freqsout = cell(length(NEWTIMEF_STRUCT.components));
%--
parfor k = 1:length(NEWTIMEF_STRUCT.components)
    %--
    fprintf('processing component %i\n',k);
    %-- temp params
    tmpp = newtimef_params;
    tmpv = NEWTIMEF_STRUCT;
    tmpd = datin;
    %--
    tmpd = tmpd(k,tmpv.pointrange,:);
    
    %-- reshape data
    timefdata  = reshape(tmpd,[1,length(tmpv.pointrange)*size(tmpd,3)]);
    
    %## RUN NEWTIMEF
    if ~isempty(timefdata)
        [P,R,mbase,ttimes,lfreqs,Pboot,Rboot,alltfX,~] ...
              = newtimef(tmpd, ...
              length(tmpv.pointrange), ...
              tmpv.timelimits, ...
              tmpv.srate, ...
              tmpv.cycles, ...
              tmpp{:});
    end

    if do_single_dat
        alltfX = single(alltfX);
    end

    datout{k}   = single(alltfX);
    timesout{k}  = ttimes;
    freqsout{k} = lfreqs;
end

%## ASSIGN DATA & ?SAVE?
ind = strcmp(EEG.subject,STUDY.datasetinfo.subject);
tinf = get_trialinfo(STUDY.datasetinfo,ind);
timef_data_struct = struct(...
    'freqs',timesout{1}, ...
    'times',freqsout{1}, ...
    'datatype','TIMEF', ...
    'datafiles',[EEG.filepath filesep EEG.filename], ...
    'datatrials',NEWTIMEF_STRUCT.indices, ...
    'parameters',newtimef_params, ...
    'trialinfo',tinf);
%-- assign data to individual fields for more memory efficient extraction
for k = 1:length(g.indices)  % for each (specified) component/channel
    timef_data_struct.(sprintf('comp%i',k)) = datout{k}; 
end

end

%% SUBFUNCTIONS ======================================================== %%
function trialinfo = get_trialinfo(datasetinfo, inds, trials)

if nargin < 1
    help std_combtrialinfo;
    return;
end

if nargin < 3
    trials = ones(1, length(datasetinfo));
end

% Inds or subjectname
if isnumeric(inds)
    if ~isempty(find(~ismember(inds(:), [datasetinfo.index]))) %#ok<EFIND>
        error(['std_combtrialinfo error: Indices '' ' num2str(inds) ' '' aout of range']);
    end
elseif ischar(inds)
    if ~ismember(inds, {datasetinfo.subject})
        error(['std_combtrialinfo error: Subject name '' ' inds ' '' is invalid']);
    else
        % Looking for indices of subject 'ind'
        inds = find(cellfun(@(x) strcmp(x,inds), {datasetinfo.subject}));
    end
    
end

if ~isfield(datasetinfo, 'trialinfo')
    trialinfo = struct([]);
    nvals = [ 1 cumsum( [trials( [ datasetinfo(inds).index ]) ])+ 1 ];
else
    % check if duration field is present
    for iDat = inds(:)'
        if ~isfield(datasetinfo(iDat).trialinfo, 'duration')
            datasetinfo(iDat).trialinfo(1).duration = [];
        end
    end
    try
        trialinfo = [ datasetinfo(inds(:)').trialinfo ];
    catch
        warning( [ 'Trial information (event) field differ for subject ' datasetinfo(inds(1)).subject '...attempting to create missing event fields' ]);
        % creating missing fields
        allFields = {};
        for iInd = inds(:)'
            allFields = union(allFields, fieldnames(datasetinfo(iInd).trialinfo));
        end
        for iInd = inds(:)'
            for iField = 1:length(allFields)
                if ~isfield(datasetinfo(iInd).trialinfo, allFields{iField})
                    datasetinfo(iInd).trialinfo(1).(allFields{iField}) = [];
                end
            end
        end
        trialinfo = [ datasetinfo(inds(:)').trialinfo ];
    end
    nvals     = [ 1 cumsum(cellfun(@length, { datasetinfo(inds).trialinfo }))+1 ];
end

fields = fieldnames(datasetinfo);
fields = setdiff( fields, { 'filepath'  'filename' 'comps' 'trialinfo' });
for iDat = 1:length(inds)
    
    for iField = 1:length(fields)
        [trialinfo(nvals(iDat):nvals(iDat+1)-1).(fields{iField})] = deal( datasetinfo(inds(iDat)).(fields{iField}) );
    end
end
end

%## =====================================================================
function std_savedat( tmpfile, structure)

    delims = find( tmpfile == '.');
    if ~isfield(structure, 'datafile') && ~isfield(structure, 'datafiles') 
        structure.datafile = [ tmpfile(1:delims(end)-1) '.set' ];
    end
    structure.date = char(datetime('today'));
    
    % fix reading problem (bug 764)
    tmpfile2  = which(tmpfile);
    if isempty(tmpfile2), tmpfile2 = tmpfile; end  
    tmpfile = tmpfile2;
    
    eeglab_options;
    if option_saveversion6
        try
            save('-v6' , tmpfile, '-struct', 'structure');
        catch
            fields = fieldnames(structure);
            for i=1:length(fields)
                eval([ fields{i} '=structure.'  fields{i} ';']);
            end
            save('-mat', tmpfile, fields{:});
        end
    else
        save('-v7.3' , tmpfile, '-struct', 'structure');
    end
end


