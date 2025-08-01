function [ersp_nolog_tf,etc_inf] = spca_get_ersp_dat(icatimef_fpath,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## Define Parser
condition_char = [];
DEF_MEAN_STRUCT = struct('do_return',true, ...
    'fcn',@(x) squeeze(mean(abs(x),3)));
% do_return_mean = true;
% return_mean_fcn = @(x) squeeze(mean(abs(x),3));
p = inputParser;
%## REQUIRED
addRequired(p,'icatimef_fpath',@ischar);
%## OPTIONAL
addOptional(p,'condition_char',condition_char,@(x) ischar(x) || isempty(x));
%## PARAMETER
addParameter(p,'MEAN_STRUCT',DEF_MEAN_STRUCT,@(x) validate_struct(x,DEF_MEAN_STRUCT));
%-- parse
parse(p,icatimef_fpath,varargin{:});
%## SET DEFAULTS
condition_char = p.Results.condition_char;
%--
MEAN_STRUCT = p.Results.MEAN_STRUCT;
MEAN_STRUCT = set_defaults_struct(MEAN_STRUCT,DEF_MEAN_STRUCT);
%% (ALTERNATIVE) USE TIMEWARPPING EEGLAB =================================== %%
% functions: newtimef.m, timefreq.m, timewarp.m
% (12/9/2023) JS, the timewarp function is the real meat that is needed to
% warp gait events to a grand average. Could use the morlet_transform_fast
% after that.
%{
%## LOAD ERSP .ICATIMEF
tf_mat_struct = load(icatimef_fpath,'-mat');

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
%}
%## GET ICATIMEF data
etc_inf = load(icatimef_fpath,'-mat', ...
    'freqs', ...
    'times', ...
    'parameters', ...
    'datatype', ...
    'datafiles', ...
    'datatrials', ...
    'trialinfo', ...
    'date');

%## NEW MORE MEMORY EFFICIENT METHOD
tt = tic();
%-- get .mat file properties
matprops = matfile(icatimef_fpath);
s = dir(icatimef_fpath);
fprintf('.icatimef data is %0.2f Gigabytes...\n',s.bytes/1e9)
%-- init extraction
fn = fieldnames(matprops);
inds = find(contains(fieldnames(matprops),'comp'));
sz = size(matprops.(fn{inds(1)}));
assd = length(sz);
%-- loops
if MEAN_STRUCT.do_return
    %-- get size
    i = 1;
    tmp = load(icatimef_fpath,'-mat',...
                        fn{inds(i)});
    val = tmp.(fn{inds(i)});
    val = MEAN_STRUCT.fcn(val);
    sz = size(val);
    assd = length(sz);
    ersp_nolog_tf = zeros([sz,length(inds)],'single');
    %-- extract data
    fprintf('Loading components (n/%i):',length(inds))
    for i = 1:length(inds)
        % fprintf(' %s',fn{inds(i)});
        fprintf(' %i',i);
        tmp = load(icatimef_fpath,'-mat',...
                        fn{inds(i)});
        val = tmp.(fn{inds(i)});
        ind = regexp(fn{inds(i)},'comp(\d*)','tokens');
        ind = double(string(ind{1}{1}));
        %-- extract specific condition
        if ~isempty(condition_char)
            cc = strcmp({etc_inf.trialinfo.cond},condition_char);
            val = squeeze(val(:,:,cc));
        end
        %-- mean fcn inline with spca protocol
        switch assd
            case 1
                ersp_nolog_tf(:,ind) = MEAN_STRUCT.fcn(val);
            case 2
                ersp_nolog_tf(:,:,ind) = MEAN_STRUCT.fcn(val);
            case 3
                ersp_nolog_tf(:,:,:,ind) = MEAN_STRUCT.fcn(val);
            case 4
                fprintf('Watcha'' doing here?\n');
                ersp_nolog_tf(:,:,:,:,ind) = MEAN_STRUCT.fcn(val);
        end
    end
    fprintf('done.\n');
else
    ersp_nolog_tf = zeros([sz,length(inds)],'complex single');
    %-- extract data
    fprintf('Loading components: (n/%i):',length(inds))
    for i = 1:length(inds)
        % fprintf(' %s.\n',fn{inds(i)})
        fprintf(' %i',i);
        val = load(icatimef_fpath,'-mat',...
                        fn{inds(i)});
        val = val.(fn{inds(i)});
        ind = regexp(fn{inds(i)},'comp(\d*)','tokens');
        ind = double(string(ind{1}{1}));
        ersp_nolog_tf(:,:,:,ind) = val;
        
        % switch assd
        %     case 1
        %         ersp_nolog_tf(:,ind) = val;
        %     case 2
        %         ersp_nolog_tf(:,:,ind) = val;
        %     case 3
        %         ersp_nolog_tf(:,:,:,ind) = val;
        %     case 4
        %         fprintf('Watcha'' doing here?\n');
        %         ersp_nolog_tf(:,:,:,:,ind) = val;
        % end
        %(03/13/20250 JS, commenting out for now as its a redundacy for
        %this extraction of data.
    end
    %-- extract specific condition
    if ~isempty(condition_char)
        cc = strcmp({tmp_etc.trialinfo.cond},condition_char);
        ersp_nolog_tf = squeeze(ersp_nolog_tf(:,:,cc,:));
    end
end
%-- time
fprintf('done. spca_get_ersp_data.m: %0.2f min\n\n',toc(tt)/(60));
end

