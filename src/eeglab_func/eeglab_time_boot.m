function [bs_tfd,bs_mask_tfd,bs_pval_tfd] = eeglab_time_boot(tf_data, ...
    varargin)
%## INPUTS
%--
% IN:
%--

tt = tic();
DEF_BOOT_STRUCT = struct(...
    'niters',1000, ...
    'alpha',0.05, ...
    'cluster_thresh',500);
%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'tf_data',@(x) isnumeric(x) && length(size(x))==3);
%## PARAMETERS
addParameter(p,'BOOT_STRUCT',DEF_BOOT_STRUCT,@(x) validate_struct(x,DEF_BOOT_STRUCT));
%## PARSE
parse(p,tf_data,varargin{:});
%## SET DEFAULTS
BOOT_STRUCT = p.Results.BOOT_STRUCT;
BOOT_STRUCT = set_defaults_struct(BOOT_STRUCT,DEF_BOOT_STRUCT);

%%
tmp_tf_mean = mean(tf_data,3);
boot_freq = 1:size(tf_data,1);
boot_time = 1:size(tf_data,2);
boot_subj = 1:size(tf_data,3);
boot_surro = zeros(size(tf_data,1),size(tf_data,2),BOOT_STRUCT.niters);

%## RANDOMIZE TIME
% scramble time samples and calculate the average across
% all times and all frequencies and store that value.
fprintf('Randomizing across time.\n')
tmp_surro = zeros(size(tf_data,1),size(tf_data,2),BOOT_STRUCT.niters);
for n = 1:BOOT_STRUCT.niters
    if mod(n,50) == 0
        fprintf('.');
    end
    %-- pick random times
    tind = randi(size(tf_data,2),[size(tf_data,2),1]);
    %-- get mean across subjs for all freqs & random times
    % tmp = mean(tf_data(boot_freq,tind,boot_subj),3);
    tmp = median(tf_data(boot_freq,tind,boot_subj),3);
    tmp_surro(:,:,n) = tmp;
end

%## RANDOMIZE SUBJECT
% Pull length(subject) surrogate averages from distribution then calc mean 
% across surrogates
fprintf('\nRandomizing across subject.\n')
for n = 1:BOOT_STRUCT.niters
    if mod(n,50) == 0
        fprintf('.');
    end
    % nind  = randi(size(tmp_surro,3),[size(tmp_surro,3),1]);
    nind  = randi(size(tmp_surro,3),[length(boot_subj),1]);
    % tmp = mean(tmp_surro(:,:,nind),3);
    tmp = median(tmp_surro(:,:,nind),3);
    boot_surro(:,:,n) = tmp;
end

%## COMPARE MEAN TF TO BOOTSTRAP
tf_pvals = stat_surrogate_pvals(boot_surro,tmp_tf_mean,'both');
% tf_pvals = stat_surrogate_pvals(boot_surro,tmp_tf_mean,'one');
%-- correct using false discovery rate
tf_pvals(tf_pvals>1)=1;
[p_masked, ~, ~, adj_p] = fdr_bh(tf_pvals,BOOT_STRUCT.alpha,'pdep',1);

%## DBSCAN IMPLEMENTATION
%{
X = zeros(length(boot_time)*length(boot_freq),2);
tmp = mean(tf_data,3);
cnt = 1;
for ti = 1:length(boot_time)
    for fi = 1:length(boot_freq)
        % X(cnt,2) = p_masked(fi,ti); %tmp(fi,ti);
        X(cnt,1) = tmp(fi,ti);
        % X(cnt,1) = adj_p(fi,ti); 
        % X(cnt,2) = log10(fi)*ti;
        X(cnt,2) = fi;
        X(cnt,3) = ti;
        cnt = cnt + 1;
    end
end
minpts = BOOT_STRUCT.cluster_thresh;
% X = squeeze(mean(tf_data,3));
% X = p_masked;
% X = adj_p;
kD = pdist2(X,X,'euc','Smallest',minpts);
%--
figure;
plot(sort(kD(end,:)));
title('k-distance graph')
xlabel('Points sorted with 50th nearest distances')
ylabel('50th nearest distances')
epsilon = 7; %0.95;
% epsilon = [];
labels = dbscan(X,epsilon,BOOT_STRUCT.cluster_thresh);
% labels = dbscan(kD,epsilon,BOOT_STRUCT.cluster_thresh);
%--
numGroups = length(unique(labels));
figure;
gscatter(X(:,1),X(:,2),labels,hsv(numGroups));
title('epsilon = 2 and minpts = 50')
grid
%}

%## BWLABELLN IMPLEMENTATION
%-- find significances that are adjacent to each other (clustered)
[lab_map,~] = bwlabeln(p_masked,2^2);
tlab_map = sort(lab_map(:),'descend');
%-- find the number of txf pts in each cluster
if ~all(unique(tlab_map) == 0)
    [occ,idx,~] = histcounts(tlab_map,unique(tlab_map));
    k_mask = ismember(lab_map,idx(occ<BOOT_STRUCT.cluster_thresh));
    fprintf('found %i unique clusters.\n',length(unique(tlab_map(k_mask))));
    k_mask = p_masked-k_mask;
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
fprintf('\ndone. eeglab_time_boot.m took %0.2f min\n',toc(tt)/60);
end
%% SUBFUNCTIONS