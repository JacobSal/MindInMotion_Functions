function [psd_log_f,etc_inf] = spca_get_psd_dat(icaspec_fpath,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tt = tic();
%## DEFAULTS
condition_char = [];
DEF_MEAN_STRUCT = struct('do_return',false, ...
    'fcn',@(x) squeeze(mean(abs(x),3)));

%## PARSER
p = inputParser;
%-- required
addRequired(p,'icaspec_fpath',@ischar);
%-- optional
addOptional(p,'condition_char',condition_char,@(x) ischar(x) || isempty(x));
%-- parameter
addParameter(p,'MEAN_STRUCT',DEF_MEAN_STRUCT,@(x) validate_struct(x,DEF_MEAN_STRUCT));
%-- parse
parse(p,icaspec_fpath,varargin{:});
%## SET DEFAULTS
condition_char = p.Results.condition_char;
%--
MEAN_STRUCT = p.Results.MEAN_STRUCT;
MEAN_STRUCT = set_defaults_struct(MEAN_STRUCT,DEF_MEAN_STRUCT);
%% (PSD) =============================================================== %%
%## GET ICATIMEF data
etc_inf = load(icaspec_fpath,'-mat', ...
    'freqs', ...
    'parameters', ...
    'datatype', ...
    'datafiles', ...
    'datatrials', ...
    'trialinfo', ...
    'date');
%-- get .mat file properties
matprops = matfile(icaspec_fpath);
s = dir(icaspec_fpath);
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
    tmp = load(icaspec_fpath,'-mat',...
                        fn{inds(i)});
    val = tmp.(fn{inds(i)});
    val = MEAN_STRUCT.fcn(val);
    sz = size(val);
    assd = length(sz);
    psd_log_f = zeros([sz,length(inds)],'single');
    %-- extract data
    fprintf('Loading components (n/%i):',length(inds))
    for i = 1:length(inds)
        % fprintf(' %s',fn{inds(i)});
        fprintf(' %i',i);
        tmp = load(icaspec_fpath,'-mat',...
                        fn{inds(i)});
        val = tmp.(fn{inds(i)});
        ind = regexp(fn{inds(i)},'comp(\d*)','tokens');
        ind = double(string(ind{1}{1}));
        %-- extract specific condition
        if ~isempty(condition_char)
            cc = strcmp({etc_inf.trialinfo.cond},condition_char);
            val = squeeze(val(:,cc));
        end
        %-- mean fcn inline with spca protocol
        switch assd
            case 1
                psd_log_f(:,ind) = MEAN_STRUCT.fcn(val);
            case 2
                psd_log_f(:,:,ind) = MEAN_STRUCT.fcn(val);
        end
    end
    fprintf('done.\n');
else
    psd_log_f = zeros([sz,length(inds)],'single');
    %-- extract data
    fprintf('Loading components: (n/%i):',length(inds))
    for i = 1:length(inds)
        % fprintf(' %s.\n',fn{inds(i)})
        fprintf(' %i',i);
        val = load(icaspec_fpath,'-mat',...
                        fn{inds(i)});
        val = val.(fn{inds(i)});
        ind = regexp(fn{inds(i)},'comp(\d*)','tokens');
        ind = double(string(ind{1}{1}));
        psd_log_f(:,:,ind) = val;

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
        cc = strcmp({etc_inf.trialinfo.cond},condition_char);
        psd_log_f = squeeze(psd_log_f(:,cc,:));
    end
end
%-- time
fprintf('\ndone. spca_get_psd_dat.m: %0.2f min\n\n',toc(tt)/(60));
end

