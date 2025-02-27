function [EEG,rej_struct_out,rmv_subj_flag] = mim_reject_subj_ics(EEG,ic_rej_in,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%- logo
cat_logo();
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
tt = tic;
%## DEFINE DEFAULTS
DEF_CHK_STRUCT = struct( ...
    'do_rmv_comps',true, ...
    'brain_thresh_int',8, ...
    'brain_exempt_int',9, ...
    'cond_chk',struct('do_reject',true, ...
        'conds',{{'0p25','0p5','0p75','1p0','rest'}}), ... % ,'flat','high','low','med'
    'ic_cut_off',5); 
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'ic_rej_in',@isstruct);
%## PARAMETER
addParameter(p,'CHK_STRUCT',DEF_CHK_STRUCT,@(x) validate_struct(x,DEF_CHK_STRUCT));
parse(p,EEG,ic_rej_in,varargin{:});
%## SET DEFAULTS
CHK_STRUCT = p.Results.CHK_STRUCT;
%-
CHK_STRUCT = set_defaults_struct(CHK_STRUCT,DEF_CHK_STRUCT);
cond_chk = CHK_STRUCT.cond_chk;
%- additional chks
if CHK_STRUCT.ic_cut_off < 2
    error('CHK_STRUCT.ic_cut_off must be greater than 1');
end
%% ===================================================================== %%
rej_struct_out = struct('subj_str',{''}, ...
    'brain_ics_count',0, ...
    'brain_ics_scores',{''}, ...
    'powpow_rej_count',0, ...
    'all_rej_count',0, ...
    'powpow_rej_nums',{''}, ...
    'brain_ic_nums',{''}, ...
    'all_rej_nums',{''}, ...
    'good_diplocs_orig',{''}, ...
    'bad_diplocs_orig',{''}, ...
    'bad_diplocs_mni',{''}, ...
    'good_diplocs_mni',{''}, ...
    'rej_notes',{'all good'});
rmv_subj_flag = false;

%## DETERMINE GOOD ICS
chk = (ic_rej_in.IC_all_brain >= CHK_STRUCT.brain_thresh_int & ...
    ic_rej_in.IC_all_brain ~= CHK_STRUCT.brain_exempt_int);
[chk_w_powpow,inds] = setdiff(find(chk),ic_rej_in.IC_powpow_rej);
tmp_bad = setdiff(find((1:size(EEG.icaweights,1))),chk_w_powpow);
tmp_good_scores = ic_rej_in.IC_all_brain(chk);
tmp_good_scores = tmp_good_scores(inds)';
tmp_good = chk_w_powpow';

%## REJECT INF STRUCT PRINTS
rej_struct_out.subj_str = EEG.subject;
%- ic criteria prints
rej_struct_out.brain_ics_count = length(tmp_good);
rej_struct_out.powpow_rej_count = length(ic_rej_in.IC_powpow_rej);
rej_struct_out.all_rej_count = length(tmp_bad);

%## REJECT INF STRUCT PRINTS
tmp = ic_rej_in.IC_powpow_rej;
if ~isempty(tmp)
    rej_struct_out.powpow_rej_nums = ['[' char(strjoin(string(tmp),',')) ']'];
end
if ~isempty(tmp_good)
    rej_struct_out.brain_ic_nums = ['[' char(strjoin(string(tmp_good),',')) ']']; 
end
if ~isempty(tmp_bad)
    rej_struct_out.all_rej_nums = ['[' char(strjoin(string(tmp_bad),',')) ']']; 
end    
if ~isempty(tmp_good_scores)
    rej_struct_out.brain_ics_scores = ['[' char(strjoin(string(tmp_good_scores),',')) ']'];
end

%## CHK CONDITIONS
if cond_chk.do_reject
    conds = unique({EEG.event.cond});
    conds = conds(~cellfun(@isempty,conds));
    chk = cellfun(@(x) any(strcmp(x,conds)),cond_chk.conds);
    if ~all(chk)
        rmv_subj_flag = true;
        str = sprintf('%s) WARNING. missing conditions. subject has only conditions: %s\n',EEG.subject,strjoin(conds,','));
        rej_struct_out.rej_notes = str;
        fprintf(str);
        return
    else
        fprintf('%s) subject has all conditions\n',EEG.subject);
    end
end
%## REMOVE & LOG THE REJECT IC PROCEDURE
EEG.etc.urreject = [];
EEG.etc.urreject.crit = [];
EEG.etc.urreject.ic_keep = [];
EEG.etc.urreject.ic_rej = [];
EEG.etc.urreject.dipfit = [];
if isempty(tmp_good)
    fprintf('%s) subject has 0 brain components\n',EEG.subject);
    rej_struct_out.rej_notes = 'no brain components post rejection';
else
    EEG.etc.urreject.crit = ic_rej_in;
    EEG.etc.urreject.ic_keep = tmp_good;
    EEG.etc.urreject.ic_rej = tmp_bad;
    EEG.etc.urreject.dipfit = EEG.dipfit;
    fprintf('%s) subject has %i brain components\n',EEG.subject, length(tmp_good));
end

%## REMOVE ICS FROM EEG
chk = length(EEG.etc.urreject.ic_keep) < CHK_STRUCT.ic_cut_off;
if chk
    fprintf('** Subject %s rejected.\n',EEG.subject);
    rej_struct_out.rej_notes = sprintf('less than %i brain components',CHK_STRUCT.ic_cut_off);
    rmv_subj_flag = true;
elseif CHK_STRUCT.do_rmv_comps
    %- (07/06/2023) JS, code taken from EEGLAB
    fprintf('%s) making dipfit & ica modifications...\n',EEG.subject)
    components = EEG.etc.urreject.ic_rej;
    fprintf('Computing projection and removing %d components ....\n', length(components));
    component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
    compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
    compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
    EEG.data(EEG.icachansind,:,:) = compproj;
    EEG.setname = [ EEG.setname ' pruned with ICA'];
    EEG.icaact  = [];
    goodinds    = setdiff_bc(1:size(EEG.icaweights,1), components);
    EEG.icawinv     = EEG.icawinv(:,goodinds);
    EEG.icaweights  = EEG.icaweights(goodinds,:);
    EEG.specicaact  = [];
    EEG.specdata    = [];
    EEG.reject      = [];
    %- iclabel mods
    fprintf('making iclabel modifications\n');
    if isfield(EEG.etc, 'ic_classification')
        if isfield(EEG.etc.ic_classification, 'ICLabel') 
            if isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
                if ~isempty(EEG.etc.ic_classification.ICLabel.classifications)
                    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodinds,:);
                end
            end
        end
    end
    try
        old_dipfit_model = EEG.dipfit.model;
        EEG.dipfit.model = EEG.dipfit.model(goodinds);
    catch e
        fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on subject %s\n',...
             'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
    end
    %- dipfit prints
    tmp = {old_dipfit_model(tmp_good).pos_old};
    if ~isempty(tmp)
        out = cellfun(@(x) ['[',sprintf('%1.0f,%1.0f,%1.0f',x),']'],tmp,'UniformOutput',false);
        out = strjoin(out,';');
        rej_struct_out.good_diplocs_orig = out;
    end
    %-
    tmp = {old_dipfit_model(tmp_bad).pos_old};
    if ~isempty(tmp)
        out = cellfun(@(x) ['[',sprintf('%1.0f,%1.0f,%1.0f',x),']'],tmp,'UniformOutput',false);
        out = strjoin(out,';');
        rej_struct_out.bad_diplocs_orig = out;
    end
    %-
    tmp = {old_dipfit_model(tmp_bad).mnipos};
    if ~isempty(tmp)
        out = cellfun(@(x) ['[',sprintf('%1.0f,%1.0f,%1.0f',x),']'],tmp,'UniformOutput',false);
        out = strjoin(out,';');
        rej_struct_out.bad_diplocs_mni = out;
    end
    %-
    tmp = {old_dipfit_model(tmp_good).mnipos};
    if ~isempty(tmp)
        out = cellfun(@(x) ['[',sprintf('%1.0f,%1.0f,%1.0f',x),']'],tmp,'UniformOutput',false);
        out = strjoin(out,';');
        rej_struct_out.good_diplocs_mni = out;
    end
else
    fprintf('%s) Not removing components...\n',EEG.subject)
end
fprintf('done. mim_reject_subj_ics.m: %0.2g s\n',toc(tt));
end

