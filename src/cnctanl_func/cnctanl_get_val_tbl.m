function [valid_tbl] = cnctanl_get_val_tbl(STUDY,cat_fext,varargin)
%FIGS_TO_EXCEL Summary of this function goes here
%   Detailed explanation goes here
% VALID_WHITE_CHARS = {'acf','ljungbox','boxpierce','limcleod'};
% VALID_IC_CHARS = {'aic','sbc','hq'};
% VALID_STAB_CHAR = 'Stability';
% VALID_PC_CHAR = 'PercentConsistency';
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct)
addRequired(p,'cat_fext',@ischar)
%## PARAMETER
parse(p,STUDY,cat_fext,varargin{:});

%%
subj_chars = {STUDY.datasetinfo.subject};
% valid_shp = [1,length(subj_chars)];
nRows = length(subj_chars);  % Or preallocate e.g., 100
[valid_tbl] = mk_tbl(nRows);

%## ADD VARS
cnt = 1;
for s_i = 1:length(subj_chars)
    try
        tmp = par_load([STUDY.datasetinfo(s_i).filepath filesep cat_fext]);
    catch e
        fprintf('%s\n',getReport(e));
        continue;
    end
    try
        srate = tmp.srate;
    catch
        srate = 0;
    end
    CAT = tmp.CAT;
    tvalid = CAT.VALIDATION;
    ticc = CAT.IC;
    nwins = length(CAT.MODEL.AR);
    %--
    valid_tbl.subj_char(cnt) = subj_chars{s_i};
    valid_tbl.cond_char(cnt) = 'all';
    valid_tbl.nchns(cnt) = CAT.nbchan;
    valid_tbl.nwins(cnt) = nwins;
    valid_tbl.srate(cnt) = srate;
    % valid_tbl.acf_w(cnt,:) = tvalid.whitestats.acf.w;
    % valid_tbl.ljungbox_w(cnt) = tvalid.whitestats.ljungbox.w;
    % valid_tbl.boxpierce_w(cnt) = tvalid.whitestats.boxpierce.w;
    % valid_tbl.limcleod_w(cnt) = tvalid.whitestats.limcleod.w;
    % valid_tbl.aic(cnt) = ticc.aic.minic; %ticc.aic.ic;
    % valid_tbl.sbc(cnt) = ticc.sbc.minic;
    % valid_tbl.hq(cnt) = ticc.hq.minic;
    % valid_tbl.aic_mo(cnt) = ticc.aic.popt; %ticc.aic.ic;
    % valid_tbl.sbc_mo(cnt) = ticc.sbc.popt;
    % valid_tbl.hq_mo(cnt) = ticc.hq.popt;
    % valid_tbl.stability(cnt) = tvalid.stabilitystats.stability;
    % valid_tbl.pctconsitency(cnt) = tvalid.PCstats.PC;
    
    %## MEAN
    valid_tbl.acf_w_mu(cnt) = mean(tvalid.whitestats.acf.w);
    valid_tbl.ljungbox_w_mu(cnt) = mean(tvalid.whitestats.ljungbox.w);
    valid_tbl.boxpierce_w_mu(cnt) = mean(tvalid.whitestats.boxpierce.w);
    valid_tbl.limcleod_w_mu(cnt) = mean(tvalid.whitestats.limcleod.w);
    valid_tbl.ljungbox_raw_mu(cnt) = mean(chk_prop(tvalid,{'whitestats','ljungbox','value'},0));
    valid_tbl.boxpierce_raw_mu(cnt) = mean(chk_prop(tvalid,{'whitestats','boxpierce','value'},0));
    tmp = chk_prop(tvalid,{'whitestats','acf','acffun'},0);
    tmp = squeeze(mean(cat(3,tmp{:}),[2,1])); % mean across channels then time
    valid_tbl.acf_raw_mu(cnt) = mean(tmp);
    valid_tbl.limcleod_raw_mu(cnt) = mean(chk_prop(tvalid,{'whitestats','limcleod','value'},0));
    valid_tbl.aic_mu(cnt) = mean(ticc.aic.minic); %ticc.aic.ic;
    valid_tbl.sbc_mu(cnt) = mean(chk_prop(ticc,{'sbc','minic'},0));
    valid_tbl.sbc_mo_mu(cnt) = mean(chk_prop(ticc,{'sbc','popt'},0));
    valid_tbl.hq_mu(cnt) = mean(ticc.hq.minic);
    valid_tbl.aic_mo_mu(cnt) = mean(ticc.aic.popt); %ticc.aic.ic;    
    valid_tbl.hq_mo_mu(cnt) = mean(ticc.hq.popt);
    %--
    tmp = chk_prop(tvalid,{'stabilitystats','lambda'},0);
    tmp = reshape(tmp,[size(tmp,1),CAT.nbchan,CAT.MODEL.morder]);
    tmp = squeeze(mean(max(tmp,[],3),2));
    tmp = squeeze(mean(tmp,2)); % mean across channels then time    
    valid_tbl.stability_mu(cnt) = mean(tmp);
    valid_tbl.pctconsitency_mu(cnt) = mean(tvalid.PCstats.PC);
    
    %## MIN
    valid_tbl.acf_w_min(cnt) = min(tvalid.whitestats.acf.w);
    valid_tbl.ljungbox_w_min(cnt) = min(tvalid.whitestats.ljungbox.w);
    valid_tbl.boxpierce_w_min(cnt) = min(tvalid.whitestats.boxpierce.w);
    valid_tbl.limcleod_w_min(cnt) = min(tvalid.whitestats.limcleod.w);
    valid_tbl.ljungbox_raw_min(cnt) = min(chk_prop(tvalid,{'whitestats','ljungbox','value'},0));
    valid_tbl.boxpierce_raw_min(cnt) = min(chk_prop(tvalid,{'whitestats','boxpierce','value'},0));
    tmp = chk_prop(tvalid,{'whitestats','acf','acffun'},0);
    tmp = squeeze(mean(cat(3,tmp{:}),[2,1])); % mean across channels then time
    valid_tbl.acf_raw_min(cnt) = min(tmp);
    valid_tbl.limcleod_raw_min(cnt) = min(chk_prop(tvalid,{'whitestats','limcleod','value'},0));
    valid_tbl.sbc_min(cnt) = min(chk_prop(ticc,{'aic','minic'},0));
    valid_tbl.sbc_min(cnt) = min(chk_prop(ticc,{'sbc','minic'},0));
    valid_tbl.sbc_mo_min(cnt) = min(chk_prop(ticc,{'sbc','popt'},0));
    valid_tbl.hq_min(cnt) = min(ticc.hq.minic);
    valid_tbl.aic_mo_min(cnt) = min(ticc.aic.popt); %ticc.aic.ic;
    valid_tbl.hq_mo_min(cnt) = min(ticc.hq.popt);
    
    %## MAX
    tmp = chk_prop(tvalid,{'stabilitystats','lambda'},0);
    tmp = reshape(tmp,[size(tmp,1),CAT.nbchan,CAT.MODEL.morder]);
    tmp = squeeze(mean(max(tmp,[],3),2));
    tmp = squeeze(mean(tmp,2)); % mean across channels then time
    valid_tbl.stability_min(cnt) = min(tmp);
    valid_tbl.pctconsitency_min(cnt) = min(tvalid.PCstats.PC);
    %--
    valid_tbl.acf_w_max(cnt) = max(tvalid.whitestats.acf.w);
    valid_tbl.ljungbox_w_max(cnt) = max(tvalid.whitestats.ljungbox.w);
    valid_tbl.boxpierce_w_max(cnt) = max(tvalid.whitestats.boxpierce.w);
    valid_tbl.limcleod_w_max(cnt) = max(tvalid.whitestats.limcleod.w);
    valid_tbl.ljungbox_raw_max(cnt) = max(chk_prop(tvalid,{'whitestats','ljungbox','value'},0));
    valid_tbl.boxpierce_raw_max(cnt) = max(chk_prop(tvalid,{'whitestats','boxpierce','value'},0));
    tmp = chk_prop(tvalid,{'whitestats','acf','acffun'},0);
    tmp = squeeze(mean(cat(3,tmp{:}),[2,1])); % mean across channels then time
    valid_tbl.acf_raw_max(cnt) = max(tmp);
    valid_tbl.limcleod_raw_max(cnt) = max(chk_prop(tvalid,{'whitestats','limcleod','value'},0));
    valid_tbl.aic_max(cnt) = max(ticc.aic.minic); %ticc.aic.ic;
    valid_tbl.sbc_max(cnt) = max(chk_prop(ticc,{'sbc','minic'},0));
    valid_tbl.sbc_mo_max(cnt) = max(chk_prop(ticc,{'sbc','popt'},0));
    valid_tbl.hq_max(cnt) = max(ticc.hq.minic);
    valid_tbl.aic_mo_max(cnt) = max(ticc.aic.popt); %ticc.aic.ic;
    valid_tbl.hq_mo_max(cnt) = max(ticc.hq.popt);
    tmp = chk_prop(tvalid,{'stabilitystats','lambda'},0);
    tmp = reshape(tmp,[size(tmp,1),CAT.nbchan,CAT.MODEL.morder]);
    tmp = squeeze(mean(max(tmp,[],3),2));
    tmp = squeeze(mean(tmp,2)); % mean across channels then time
    valid_tbl.stability_max(cnt) = max(tmp);
    valid_tbl.pctconsitency_max(cnt) = max(tvalid.PCstats.PC);
    
    %## STD
    valid_tbl.acf_w_std(cnt) = std(tvalid.whitestats.acf.w);
    valid_tbl.ljungbox_w_std(cnt) = std(tvalid.whitestats.ljungbox.w);
    valid_tbl.boxpierce_w_std(cnt) = std(tvalid.whitestats.boxpierce.w);
    valid_tbl.limcleod_w_std(cnt) = std(tvalid.whitestats.limcleod.w);
    valid_tbl.ljungbox_raw_std(cnt) = std(chk_prop(tvalid,{'whitestats','ljungbox','value'},0));
    valid_tbl.boxpierce_raw_std(cnt) = std(chk_prop(tvalid,{'whitestats','boxpierce','value'},0));
    tmp = chk_prop(tvalid,{'whitestats','acf','acffun'},0);
    tmp = squeeze(mean(cat(3,tmp{:}),[2,1])); % mean across channels then time
    valid_tbl.acf_raw_std(cnt) = std(tmp);
    valid_tbl.limcleod_raw_std(cnt) = std(chk_prop(tvalid,{'whitestats','limcleod','value'},0));
    valid_tbl.aic_std(cnt) = std(ticc.aic.minic); %ticc.aic.ic;
    valid_tbl.sbc_std(cnt) = std(chk_prop(ticc,{'sbc','minic'},0));
    valid_tbl.sbc_mo_std(cnt) = std(chk_prop(ticc,{'sbc','popt'},0));
    valid_tbl.hq_std(cnt) = std(ticc.hq.minic);
    valid_tbl.aic_mo_std(cnt) = std(ticc.aic.popt); %ticc.aic.ic;
    valid_tbl.hq_mo_std(cnt) = std(ticc.hq.popt);
    tmp = chk_prop(tvalid,{'stabilitystats','lambda'},0);
    tmp = reshape(tmp,[size(tmp,1),CAT.nbchan,CAT.MODEL.morder]);
    tmp = squeeze(mean(max(tmp,[],3),2));
    tmp = squeeze(mean(tmp,2)); % mean across channels then time
    valid_tbl.stability_std(cnt) = std(tmp);
    valid_tbl.pctconsitency_std(cnt) = std(tvalid.PCstats.PC);
    %--
    cnt = cnt + 1;
end
valid_tbl = valid_tbl(valid_tbl.nwins ~= 0,:);
%--
end
%% ===================================================================== %%
function [out_val] = chk_prop(struct,prop_names,def_val)
    try
        for p = 1:length(prop_names)
            out_val = struct.(prop_names{p});
            struct = out_val;
        end
    catch
        fprintf('Not a property: %s\n',prop_names{p});
        out_val = def_val;
    end
end

%## =======================================================================
function [tbl] = mk_tbl(nRows)
    tbl = table( ...
        repmat("", nRows, 1), ... % string
        repmat("", nRows, 1), ... % string
        zeros( nRows, 1), ...      % numeric
        zeros( nRows, 1), ...       % numeric
        zeros( nRows, 1), ...       %--
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...    
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...     
        zeros( nRows, 1), ...       %--
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...    
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...     
        zeros( nRows, 1), ...       %--
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...    
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...    
        zeros( nRows, 1), ...       %--
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...    
        zeros( nRows, 1), ...       
        zeros( nRows, 1), ...     
        'VariableNames', {'subj_char', 'cond_char', 'nchns', 'nwins' ...
        'acf_w_mu','ljungbox_w_mu','boxpierce_w_mu','limcleod_w_mu', ...
        'aic_mu','sbc_mu','hq_mu','aic_mo_mu','sbc_mo_mu','hq_mo_mu', ...
        'stability_mu','pctconsitency_mu', ...
        'acf_w_std','ljungbox_w_std','boxpierce_w_std','limcleod_w_std', ...
        'aic_std','sbc_std','hq_std','aic_mo_std','sbc_mo_std','hq_mo_std', ...
        'stability_std','pctconsitency_std', ...
        'acf_w_min','ljungbox_w_min','boxpierce_w_min','limcleod_w_min', ...
        'aic_min','sbc_min','hq_min','aic_mo_min','sbc_mo_min','hq_mo_min', ...
        'stability_min','pctconsitency_min' ...
        'acf_w_max','ljungbox_w_max','boxpierce_w_max','limcleod_w_max', ...
        'aic_max','sbc_max','hq_max','aic_mo_max','sbc_mo_max','hq_mo_max', ...
        'stability_max','pctconsitency_max'} ...
    );
end
