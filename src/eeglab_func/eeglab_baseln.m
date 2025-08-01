function [tf_data_full,tf_data_crop,baseln_com,baseln_trial] = eeglab_baseln(tf_data,times,freqs,base_tlims,base_flims,...
    varargin)
%## INPUTS
%--
% IN:
%       tf_data, CELL
%           should be log-transformed time-frequency data (freq x time x subjs)
%       times, DOUBLE/SINGLE          
%           vector of times in milliseconds with size(tf_data{1},2)
%       freqs, DOUBLE/SINGLE
%           vector of frequencies in hertz with size(tf_data{1},1)
%       basetimes, DOUBLE/SINGLE
%           a pair of times(ms) [BEGINNING,END] to baseline tf_data to
%       basefreqs, DOUBLE/SINGLE
%           a pair of freqs(hz) [BEGINNING,END] to baseline tf_data to.
%           This is more for plotting and doesn't have a major effect on
%           baselining outcome.
%--
DO_COMMON_BASE = false;
DO_SUBJ_BASE = true;
GROUPWISE_PARAMS = struct('single_grp_treatment',false,...
    'group_wts',[1,1,1]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'tf_data',@iscell); % this should be shaped [freqs, times, subjs]
addRequired(p,'times',@isnumeric);
addRequired(p,'freqs',@isnumeric);
addRequired(p,'base_tlims',@isnumeric);
addRequired(p,'base_flims',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'DO_COMMON_BASE',DO_COMMON_BASE,@islogical);
addParameter(p,'DO_SUBJ_BASE',DO_SUBJ_BASE,@islogical);
addParameter(p,'GROUPWISE_PARAMS',GROUPWISE_PARAMS,@struct);
parse(p,tf_data,times,freqs,base_tlims,base_flims,varargin{:});
%## SET DEFAULTS
DO_COMMON_BASE = p.Results.DO_COMMON_BASE;
DO_SUBJ_BASE = p.Results.DO_SUBJ_BASE;
GROUPWISE_PARAMS = p.Results.GROUPWISE_PARAMS;
%%
baseln_com = zeros(size(tf_data{1},1),size(tf_data{1},3));
baseln_trial = cell(size(tf_data));
base_tinds = (times>=base_tlims(1) & times<=base_tlims(end));
base_finds = (freqs>=base_flims(1) & freqs<=base_flims(end));
%--
sz = size(tf_data{1});
indt = find(sz == length(times));
indf = find(sz == length(freqs));
inds = find(indt ~= (1:length(sz)) & indf ~= (1:length(sz)));
for i = 1:size(tf_data,1)
    for j = 1:size(tf_data,2)
        tf_data{i,j} = permute(tf_data{i,j},[indf,indt,inds]); % reshape so its in the form [freq, time, subj]
    end
end
%--
% sz = size(tf_data{1});
% indt = find(sz == length(times));
% indf = find(sz == length(freqs));
% inds = find(indt ~= (1:length(sz)) & indf ~= (1:length(sz)));
% tf_data_bcom1 = cell(size(tf_data));
% tf_data_bcom2 = cell(size(tf_data));
% tf_data_log_subBase = cell(size(tf_data));
% tf_data_btrial = cell(size(tf_data));
tmp_tf_data = tf_data;
%## PER SUBJECT, GROUP, CONDITION BASELINE
if DO_SUBJ_BASE
    %## OPT1
    tmp_tf_data = cellfun(@(x) sub_base(x,base_tinds,2),tmp_tf_data,'UniformOutput',false);
    % tmp_tf_data = cellfun(@(x) sub_base(x,base_tinds,1),tmp_tf_data,'UniformOutput',false);
    
    %## OPT2
    % for i = 1:size(tmp_tf_data,2) % groups
    %     for j = 1:size(tmp_tf_data,1) % conditions
    %         ted = tmp_tf_data{j,i}(:,:,:);
    %         %-- mean power across time for each person (or trial)
    %         tmp_base_subj = mean(ted(:,base_tinds,:),2);
    %         %-- subtract out each subjects baseline
    %         ssz = size(ted,3);
    %         for si = 1:ssz
    %             tmp = squeeze(ted(:,:,si));
    %             tmp = bsxfun(@minus,tmp,tmp_base_subj(:,si));
    %             tf_data_btrial{j,i}(:,:,si) = permute(tmp,[1,2]);
    %         end
    %         %-- store base
    %         baseln_trial{j,i} = tmp_base_subj;
    %     end
    % end
    % tmp_tf_data = tf_data_btrial;
end

%## COMMON BASELINE ACROSS GROUP
if DO_COMMON_BASE
    %--
    for i = 1:size(tmp_tf_data,2)
        tmp_log = cat(3,tmp_tf_data{:,i});
        tmp_logb = mean(tmp_log,3);
        %## OPT1
        tmp = cellfun(@(x) com_base(x,base_tinds,tmp_logb),tmp_tf_data(:,i),'UniformOutput',false);
        tmp_tf_data(:,i) = tmp;
    end

    %## OPT 2?
    % %-- loop over groups
    % for i = 1:size(tmp_tf_data,2)
    %     %-- loop over subjects
    %     for si = 1:size(tmp_tf_data{1,i},3)
    %         %-- convert log-spectrum to non-log for all conditions
    %         for j = 1:size(tmp_tf_data,1)
    %             tmp_nolog(:,:,j) = 10.^(tmp_tf_data{j,i}(:,:,si)/20);
    %         end
    %         %-- mean across conditions then mean across times
    %         tmp_base_com = mean(mean(tmp_nolog(:,base_tinds,:),3),2);
    %         %-- log-subtraction of baseline for a particular subject
    %         for j = 1:size(tmp_tf_data,1)
    %             tf_data_bcom1{j,i}(:,:,si) = tmp_nolog(:,:,j)./(repmat(tmp_base_com,[1,size(tmp_nolog,2)]));
    %             tf_data_bcom2{j,i}(:,:,si) = 20*log10(tf_data_bcom1{j,i}(:,:,si)); 
    %         end
    %         %-- store
    %         baseln_com(:,si) = tmp_base_com;
    %     end
    % end
    %-- store
    % tmp_tf_data = tf_data_bcom2;
end
tf_data_full = tmp_tf_data;
tf_data_crop = cellfun(@(x) x(base_finds,base_tinds,:),tmp_tf_data,'uniformoutput',false);
end
%% SUBFUNCTIONS
%## SUBJECT BASELINE FUNCTION
function dato = sub_base(datin,bt_inds,tdim)
    %-- mean power across time for each person (or trial)
    switch tdim
        case 1
            tmp_bs = mean(datin(bt_inds,:,:),1);
        case 2
            tmp_bs = mean(datin(:,bt_inds,:),2);
        case 3
            tmp_bs = mean(datin(:,:,bt_inds),3);
        otherwise
            error('bad time dimension.');
    end
    %-- subtract out each subjects baseline
    dato = bsxfun(@minus,datin,tmp_bs);
end

%## COMMON BASELINE FUNCTION
function dato = com_base(datin,bt_inds,bt_dat)
    % %-- mean power across time for each person (or trial)
    % switch tdim
    %     case 1
    %         tmp_bs = mean(bt_dat(bt_inds,:,:),1);
    %         tmp_bs = repmat(tmp_bs,[])
    %     case 2
    %         tmp_bs = mean(bt_dat(:,bt_inds,:),2);
    %     case 3
    %         tmp_bs = mean(bt_dat(:,:,bt_inds),3);
    %     otherwise
    %         error('bad time dimension.');
    % end
    %--
    tmp_bs = mean(bt_dat(:,bt_inds,:),2);
    tmp_bs = repmat(tmp_bs,[1,1,size(datin,3)]);
    %-- subtract out baseline data
    dato = bsxfun(@minus,datin,tmp_bs);
end