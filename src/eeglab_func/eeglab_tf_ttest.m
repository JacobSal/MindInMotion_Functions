function [bs_tfd,bs_mask_tfd,bs_pval_tfd] = eeglab_tf_ttest(tf_data, comp_data, ...
    varargin)
%## INPUTS
%--
% IN:
%--
% comp_data = 0;
tt = tic();
DEF_TTEST_STRUCT = struct(...
    'tail','both', ...
    'alpha',0.05, ...
    'cluster_thresh',300);
%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'tf_data',@(x) isnumeric(x) && length(size(x))==3);
addRequired(p,'comp_data',@(x) isnumeric(x));

%## PARAMETERS
addParameter(p,'TTEST_STRUCT',DEF_TTEST_STRUCT,@(x) validate_struct(x,DEF_TTEST_STRUCT));
%## PARSE
parse(p,tf_data,comp_data,varargin{:});
%## SET DEFAULTS
TTEST_STRUCT = p.Results.TTEST_STRUCT;
TTEST_STRUCT = set_defaults_struct(TTEST_STRUCT,DEF_TTEST_STRUCT);

%%
tmp_tf_mean = mean(tf_data,3);
tffreq = 1:size(tf_data,1);
tftime = 1:size(tf_data,2);
% tfsubj = 1:size(tf_data,3);
tf_pvals = zeros(size(tf_data,1),size(tf_data,2),1);
tstat = zeros(size(tf_data,1),size(tf_data,2),1);
tmask = zeros(size(tf_data,1),size(tf_data,2),1);

%## TTEST
%-- test if a specific time frequency bin is greater or lesser than 0
fprintf('Performing t-tests.\n');
for ti = 1:length(tftime)
    for fi = 1:length(tffreq)
        if mod(ti+fi,75)==0
            fprintf('.');
        end
        %-- test if data > 0 (right tailed?)
        % tail::: 'both' | 'right' | 'left'
        [h,p,ci,stats] = ttest(squeeze(tf_data(fi,ti,:)),comp_data, ...
            'Alpha',TTEST_STRUCT.alpha, ...
            'Tail',TTEST_STRUCT.tail); %#ok<ASGLU>
        %--
        % [h,p,ci,stats] = ttest(squeeze(tf_data(fi,ti,:)),0, ...
        %     'Alpha',0.05, ...
        %     'Tail','right'); %#ok<ASGLU>
        %--
        tf_pvals(fi,ti) = p;
        tstat(fi,ti) = stats.tstat;
        tmask(fi,ti) = h;
    end
end
fprintf('\n');

%## COMPARE MEAN TF TO BOOTSTRAP
%-- correct using false discovery rate
tf_pvals(tf_pvals>1)=1;
[p_masked, ~, ~, adj_p] = fdr_bh(tf_pvals,TTEST_STRUCT.alpha,'pdep',1);
%-- find significances that are adjacent to each other (clustered)
[lab_map,~] = bwlabeln(p_masked);
tlab_map = sort(lab_map(:),'descend');
%-- find the number of txf pts in each cluster
if ~all(unique(tlab_map) == 0)    
    [occ,idx,~] = histcounts(tlab_map,unique(tlab_map));
    k_mask = ismember(lab_map,idx(occ<TTEST_STRUCT.cluster_thresh));
    k_mask = p_masked-k_mask;
    fprintf('found %i unique clusters.\n',length(unique(tlab_map(~k_mask))));
    %-- mask mean data
    tmp = tmp_tf_mean; 
    tmp(~k_mask) = 0;
else
    tmp = zeros(size(tmp_tf_mean));
end
%-- store
bs_tfd = tmp_tf_mean;
bs_mask_tfd = tmp;
bs_pval_tfd = adj_p;
%--
fprintf('\ndone. eeglab_tf_ttest.m took %0.2f min\n',toc(tt)/60);
end
%% SUBFUNCTIONS