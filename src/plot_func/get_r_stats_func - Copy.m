function [STATS_STRUCT,CONFINT_STRUCT,ranef,raw_values,tmp_stats] = get_r_stats_func(rstats_table,cl_n,varargin)
%GET_R_STATS_FUNC Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

% params.str_offset = [-0.1,-0.05];
% params.str_font_size = 7;
% params.ci_bar_width = 0.15;
% params.ci_bar_xpos = 0; % 0.5;
% params.ci_bar_linespecs = {'LineStyle','-','LineWidth',2,'Color','k'}; %[.9 .9 .9]
% params.coeff_desm = [[1,0]; ...
%     [0,1]; ...
%     [-1,-1]];

%## TIME
tic
%## DEFINE DEFAULTS
% COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_charH2000''s','group_charH3000''s', ...
%     'speed_cond_num:group_charH2000''s','speed_cond_num:group_charH3000''s'};
% ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
% ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
% COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_charH2000''s','group_charH3000''s'};
%--
COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_char1','group_char2', ...
    'speed_cond_num:group_char1','speed_cond_num:group_char2'};
ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char1','group_char2'};
%--
DEF_EXTRACT_STRUCT = struct( ...
    'group_chars', {{'H1000''s','H2000''s','H3000''s'}}, ...
    'group_order', categorical({'1','2','3'}), ...
    'model_char_int', 'speed_group_intact_all', ...
    'model_char_group', 'speed_group_all', ...
    'group_char', 'all', ...
    'anv_chars_int', {ANV_CHARS_INT}, ...
    'anv_chars_group', {ANV_CHARS_GROUP}, ...
    'coeff_chars_int', {COEFF_CHARS_INT}, ...
    'coeff_chars_group', {COEFF_CHARS_GROUP}, ...
    'eeg_measure', {''}, ...
    'kin_measure', {'none'}, ...
    'coeff_study', {'speed'}, ...
    'str_option',{'simple'}, ...
    'str_offset', [-0.1, -0.05], ...
    'str_font_size', 7, ...
    'do_extract_raw',false, ...
    'ci_bar_width', 0.15, ...
    'ci_bar_xpos', 0, ...
    'ci_bar_linespecs', {{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
    'coeff_desm', [ ...
        1, 0; ...
        0, 1; ...
       -1,-1 ...
    ] ...
);
%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'rstats_table',@istable);
addRequired(p,'cl_n',@isnumeric);
%## PARAMETER
addParameter(p,'EXTRACT_STRUCT',DEF_EXTRACT_STRUCT,@(x) validate_struct(x,DEF_EXTRACT_STRUCT));
%##
parse(p,rstats_table,cl_n,varargin{:});
%## SET DEFAULTS
%-- extract vars
EXTRACT_STRUCT = p.Results.EXTRACT_STRUCT;
%-- set defaults
EXTRACT_STRUCT = set_defaults_struct(EXTRACT_STRUCT,DEF_EXTRACT_STRUCT);
%--
params = EXTRACT_STRUCT;

%## INITIATE OUTPUTS
ranef = [];
tssc = struct('var_type','continuous', ...
    'anova',1, ...
    'multc_pvals',[],...
    'multc_pvals_pairs',[],...
    'regress_line',[], ... 
    'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
    'regress_xvals',[],... % continuous predictors
    'order',{{}});
tssg = struct('var_type','categorical', ...
    'anova',1, ...
    'multc_pvals',[],...
    'multc_pvals_pairs',[],...
    'regress_line',[],... 
    'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
    'regress_xvals',0,... % continuous predictors
    'order',{categorical({''})});
STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
    'cond_stats',tssc, ...
    'group_stats',tssg, ...
    'stats_char',struct('str',{{''}}, ...
        'offsets',[], ...
        'font_size',[], ...
        'do_display',true));
CONFINT_STRUCT = struct('do_display',false, ...
    'y_bnds',[], ...
    'x_vals',[], ...
    'errbar_struct',struct('line_specs',{params.ci_bar_linespecs}, ...
        'err_bar_width',params.ci_bar_width));
%--
raw_values = struct('anova_p',[], ...
    'anova_s',[], ...
    'coeffs',[], ...
    'fsq',[], ...
    'rsq',[]);
%% ===================================================================== %%
%## EXTRACT STATS INFO
if ~isfield(params,'kin_measure')
    params.kin_measure = [];
end
if isempty(params.kin_measure)
    tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char_int) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
else
    tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char_int) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.kinematic_char,params.kin_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
end
tmp_stats = rstats_table(tmp_stats,:);
% tmp_stats = tmp_stats(2,:);
%--
tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
%-- anova p-values
anv_p = zeros(length(params.anv_chars_int),1);
for cc = 1:length(params.anv_chars_int)
    ind = strcmp(params.anv_chars_int{cc},tmp_ac);
    if ~isempty(ind)
        anv_p(cc) = tmp_anv(ind);              
    else
        fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
    end            
end

tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
tmp_lme_cc = strsplit(tmp_stats.lme_cc{1},',');
tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
% tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
%-- rndm inter
tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
ranef = struct('char',tmp_ranef_char, ...
    'int',tmp_ranef_n);    

%## LME PARAM EXTRACT
lme_est = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_est');
lme_stat = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_stat');
lme_pv = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_pval');
%--
anv_p = tbl2num(params.anv_chars_int,tmp_ac,tmp_stats,'anv_pvals');
anv_stat = tbl2num(params.anv_chars_int,tmp_ac,tmp_stats,'anv_stats');
coeffs = tbl2num(params.coeff_chars_int,tmp_cc,tmp_stats,'coeffs');


%## CHECK FOR SIGNIFICANT INTERACTION
if length(anv_p) > 3 && ~EXTRACT_STRUCT.do_extract_raw
    tt = {'grp','int'};
    anl_type = double(anv_p(4) > 0.05)+1;
    anl_type = tt{anl_type};
end
switch anl_type
    case 'grp'
        %## USE INTERACTION MODEL
    tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
    tmp_lme_cc = strsplit(tmp_stats.lme_cc{1},',');
    tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
    tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
    tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
    % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
    tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
    tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
    %-- rndm inter
    tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
    tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    ranef = struct('char',tmp_ranef_char, ...
        'int',tmp_ranef_n);    
    
    %## LME PARAM EXTRACT
    lme_est = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_est');
    lme_stat = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_stat');
    lme_pv = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_pval');
    %--
    anv_stat = tbl2num(params.anv_chars_int,tmp_ac,tmp_stats,'anv_stats');
    coeffs = tbl2num(params.coeff_chars_int,tmp_cc,tmp_stats,'coeffs');
    %## EXTRACT MODEL PARAMS
    %-- cohens f^2 values
    fsq_chars = strcmp(params.anv_chars_int,'(Intercept)');
    fsq_chars = params.anv_chars_int(~fsq_chars);
    fsqs = zeros(length(fsq_chars),1);
    for cc = 1:length(fsq_chars)
        ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
        if ~isempty(ind)
            fsqs(cc) = tmp_fsq(ind);
        else
            fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
        end
    end
    %-- confidence intervals
    if ~all(cellfun(@isempty,tmp_ci_chars))
        ci_chars = strcmp(params.group_chars,'(Intercept)');
        ci_chars = params.group_chars(~ci_chars);
        cis = zeros(length(ci_chars),1,2);
        for cc = 1:length(ci_chars)
            ind = strcmp(ci_chars{cc},tmp_ci_chars);
            if ~isempty(ind)
                cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
            else
                fprintf("Coefficient %s not found.\n",params.group_chars{cc})
            end
        end
        %--
        if anv_p(4) > 0.05 && anv_p(3) > 0.05
            chkd = false;
        else
            chkd = true;
        end
        CONFINT_STRUCT = struct('do_display',chkd, ...
                'y_bnds',cis, ...
                'x_vals',repmat(params.ci_bar_xpos,[3,1]), ...
                'errbar_struct',struct('line_specs',{params.ci_bar_linespecs}, ...
                    'err_bar_width',params.ci_bar_width));
    end

    
    
    switch params.coeff_study
        case 'sbs'
            %-- (04/17/2025) JS, coeff organization for step-by-step study
            regl = [[coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)+params.coeff_desm(1,1)*coeffs(5)+params.coeff_desm(1,2)*coeffs(6)]; ... %lvl 1
                [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)+params.coeff_desm(2,1)*coeffs(5)+params.coeff_desm(2,2)*coeffs(6)]; ...
                [coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)+params.coeff_desm(3,1)*coeffs(5)+params.coeff_desm(3,2)*coeffs(6)]]; %lvl 1=3
        case 'speed'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            regl = [[coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)+params.coeff_desm(3,1)*coeffs(5)+params.coeff_desm(3,2)*coeffs(6)]; ...
                [coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)+params.coeff_desm(1,1)*coeffs(5)+params.coeff_desm(1,2)*coeffs(6)]; ... %lvl 1
                [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)+params.coeff_desm(2,1)*coeffs(5)+params.coeff_desm(2,2)*coeffs(6)]]; %lvl 1=3
        case 'speed:kin'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            % regl = [[coeffs(1),coeffs(2)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4)]; ...
            %     [coeffs(1),coeffs(2)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4)]; ... %lvl 1
            %     [coeffs(1),coeffs(2)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4)]]; %lvl 1=3
            regl = [coeffs(1),coeffs(2),coeffs(3),coeffs(4)];
        case 'speed:kin:grp'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            % regl = [[coeffs(1),coeffs(2)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4)]; ...
            %     [coeffs(1),coeffs(2)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4)]; ... %lvl 1
            %     [coeffs(1),coeffs(2)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4)]]; %lvl 1=3
            regl = coeffs;
    end

    %## ASSIGN STATS
    str = {'','',''};
    switch EXTRACT_STRUCT.str_option
        case 'simple'
            [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
        case 'information'
            [str] = get_inform_str(anv_p,tmp_stats.r2_c_int,regl);
    end
    % if EXTRACT_STRUCT.do_display_str
    %     % [str] = get_anv_str(anvs,tmp_stats.r2_c_int,fsqs);
    %     [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
    %     [str] = get_inform_str(anv_p,r2s,regl);
    % end

    %## RAW OUTPUT
    %--
    raw_values = struct('anova_p',anv_p, ...
        'anova_s',anv_stat, ...
        'coeffs',coeffs, ...
        'lme_ccs',{tmp_lme_cc}, ...
        'lme_pv',lme_pv, ...
        'lme_est',lme_est, ...
        'lme_stat',lme_stat, ...
        'fsq',fsqs, ...
        'rsq',[tmp_stats.r2_c_int,tmp_stats.r2_m_int]);

    %##  
    tssc = struct('var_type','continuous', ...
        'anova',anv_p(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',regl, ...
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',(0:5)*0.25,... % continuous predictors
        'order',{{}});
    tssg = struct('var_type','categorical', ...
        'anova',anv_p(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',[],... 
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',0,... % continuous predictors
        'order',{params.group_order});
    STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
        'cond_stats',tssc, ...
        'group_stats',tssg, ...
        'stats_char',struct('str',{str}, ...
            'offsets',params.str_offset, ...
            'font_size',params.str_font_size, ...
            'do_display',true));
    case 'int'
    case 'raw'
end
if anv_p(4) > 0.05 && ~EXTRACT_STRUCT.do_extract_raw
    %## USE NON-INTERACTION MODEL
    if isempty(params.kin_measure)
        tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
            strcmp(rstats_table.model_char,params.model_char_group) &...
            strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
            strcmp(rstats_table.group_char,params.group_char);
    else
        tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
            strcmp(rstats_table.model_char,params.model_char_group) &...
            strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
            strcmp(rstats_table.kinematic_char,params.kin_measure) &...
            strcmp(rstats_table.group_char,params.group_char);
    end
    tmp_stats = rstats_table(tmp_stats,:);
    if ~isempty(tmp_stats)
        %--
        tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
        tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
        %-- anova p-values
        anv_p = zeros(length(params.anv_chars_group),1);
        for cc = 1:length(params.anv_chars_group)
            ind = strcmp(params.anv_chars_group{cc},tmp_ac);
            if ~isempty(ind)
                anv_p(cc) = tmp_anv(ind);              
            else
                fprintf("Coefficient %s not found.\n",params.anv_chars_group{cc})
            end            
        end       
        tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
        tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
        tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
        tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
        tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
        % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
        tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
        tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
        %-- rndm inter
        tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
        tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
        ranef = struct('char',tmp_ranef_char, ...
            'int',tmp_ranef_n);
        %-- model coefficients
        coeffs = zeros(length(params.coeff_chars_group),1);        
        for cc = 1:length(params.coeff_chars_group)
            ind = strcmp(params.coeff_chars_group{cc},tmp_cc);
            if ~isempty(ind)
                coeffs(cc) = tmp_coeffs(ind);              
            else
                fprintf("Coefficient %s not found.\n",params.coeff_chars_group{cc})
            end            
        end        
        %-- cohens f^2 values
        fsq_chars = strcmp(params.anv_chars_group,'(Intercept)');
        fsq_chars = params.anv_chars_group(~fsq_chars);
        fsqs = zeros(length(fsq_chars),1);
        for cc = 1:length(fsq_chars)
            ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
            if ~isempty(ind)
                fsqs(cc) = tmp_fsq(ind);
            else
                fprintf("Coefficient %s not found.\n",params.anv_chars_group{cc})
            end
        end
        %-- confidence intervals
        ci_chars = strcmp(params.group_chars,'(Intercept)');
        ci_chars = params.group_chars(~ci_chars);
        cis = zeros(length(ci_chars),1,2);
        for cc = 1:length(ci_chars)
            ind = strcmp(ci_chars{cc},tmp_ci_chars);
            if ~isempty(ind)
                cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
            else
                fprintf("Coefficient %s not found.\n",params.group_chars{cc})
            end
        end
    
        %## ASSIGN STATS
        % str = {'','',''};
        % if EXTRACT_STRUCT.do_display_str
        %     % [str] = get_anv_str(anvs,tmp_stats.r2_c_int,fsqs);
        %     [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
        % end
        %--
        if anv_p(3) > 0.05 && anv_p(2) > 0.05
            chkd = false;
        else
            chkd = true;
        end
        CONFINT_STRUCT = struct('do_display',chkd, ...
                'y_bnds',cis, ...
                'x_vals',repmat(params.ci_bar_xpos,[3,1]), ...
                'errbar_struct',struct('line_specs',{params.ci_bar_linespecs}, ...
                    'err_bar_width',params.ci_bar_width));
        
        switch params.coeff_study
            case 'sbs'
                %-- (04/17/2025) JS, coeff organization for step-by-step study
                regl = [[coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)]; ... %lvl 1
                    [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)]; ...
                    [coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)]]; %lvl 1=3
            case 'speed'
                %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
                regl = [[coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)]; ...
                    [coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)]; ... %lvl 1
                    [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)]]; %lvl 1=3
        end

        %## ASSIGN STATS
        str = {'','',''};
        switch EXTRACT_STRUCT.str_option
            case 'simple'
                [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
            case 'information'
                [str] = get_inform_str(anv_p,tmp_stats.r2_c_int,regl);
        end
        %--
        tssc = struct('var_type','continuous', ...
            'anova',anv_p(2), ...
            'multc_pvals',[],...
            'multc_pvals_pairs',[],...
            'regress_line',regl, ... 
            'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
            'regress_xvals',(0:5)*0.25,... % continuous predictors
            'order',{{}});
        tssg = struct('var_type','categorical', ...
            'anova',anv_p(3), ...
            'multc_pvals',[],...
            'multc_pvals_pairs',[],...
            'regress_line',[],... 
            'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
            'regress_xvals',0,... % continuous predictors
            'order',{params.group_order});
        STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
            'cond_stats',tssc, ...
            'group_stats',tssg, ...
            'stats_char',struct('str',{str}, ...
                'offsets',params.str_offset, ...
                'font_size',params.str_font_size, ...
                'do_display',true));
    else
        fprintf('No non-interaction stats found... Returning empty structs\n');        
    end
elseif anv_p(4) < 0.05 || EXTRACT_STRUCT.do_extract_raw
    %## USE INTERACTION MODEL
    tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
    tmp_lme_cc = strsplit(tmp_stats.lme_cc{1},',');
    tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
    tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
    tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
    % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
    tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
    tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
    %-- rndm inter
    tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
    tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    ranef = struct('char',tmp_ranef_char, ...
        'int',tmp_ranef_n);    
    
    %## LME PARAM EXTRACT
    lme_est = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_est');
    lme_stat = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_stat');
    lme_pv = tbl2num(params.coeff_chars_int,tmp_lme_cc,tmp_stats,'lme_pval');
    %--
    anv_stat = tbl2num(params.anv_chars_int,tmp_ac,tmp_stats,'anv_stats');
    coeffs = tbl2num(params.coeff_chars_int,tmp_cc,tmp_stats,'coeffs');
    %## EXTRACT MODEL PARAMS
    %-- anova STATS
    % anv_stat = zeros(length(params.anv_chars_int),1);
    % for cc = 1:length(params.anv_chars_int)
    %     ind = strcmp(params.anv_chars_int{cc},tmp_ac);
    %     if ~isempty(ind)
    %         anv_stat(cc) = tmp_anv_s(ind);              
    %     else
    %         fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
    %     end            
    % end
    %-- model coefficients
    % coeffs = zeros(length(params.coeff_chars_int),1);        
    % for cc = 1:length(params.coeff_chars_int)
    %     ind = strcmp(params.coeff_chars_int{cc},tmp_cc);
    %     if ~isempty(ind)
    %         coeffs(cc) = tmp_coeffs(ind);              
    %     else
    %         fprintf("Coefficient %s not found.\n",params.coeff_chars_group{cc})
    %     end            
    % end  
    %-- cohens f^2 values
    fsq_chars = strcmp(params.anv_chars_int,'(Intercept)');
    fsq_chars = params.anv_chars_int(~fsq_chars);
    fsqs = zeros(length(fsq_chars),1);
    for cc = 1:length(fsq_chars)
        ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
        if ~isempty(ind)
            fsqs(cc) = tmp_fsq(ind);
        else
            fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
        end
    end
    %-- confidence intervals
    if ~all(cellfun(@isempty,tmp_ci_chars))
        ci_chars = strcmp(params.group_chars,'(Intercept)');
        ci_chars = params.group_chars(~ci_chars);
        cis = zeros(length(ci_chars),1,2);
        for cc = 1:length(ci_chars)
            ind = strcmp(ci_chars{cc},tmp_ci_chars);
            if ~isempty(ind)
                cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
            else
                fprintf("Coefficient %s not found.\n",params.group_chars{cc})
            end
        end
        %--
        if anv_p(4) > 0.05 && anv_p(3) > 0.05
            chkd = false;
        else
            chkd = true;
        end
        CONFINT_STRUCT = struct('do_display',chkd, ...
                'y_bnds',cis, ...
                'x_vals',repmat(params.ci_bar_xpos,[3,1]), ...
                'errbar_struct',struct('line_specs',{params.ci_bar_linespecs}, ...
                    'err_bar_width',params.ci_bar_width));
    end

    
    
    switch params.coeff_study
        case 'sbs'
            %-- (04/17/2025) JS, coeff organization for step-by-step study
            regl = [[coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)+params.coeff_desm(1,1)*coeffs(5)+params.coeff_desm(1,2)*coeffs(6)]; ... %lvl 1
                [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)+params.coeff_desm(2,1)*coeffs(5)+params.coeff_desm(2,2)*coeffs(6)]; ...
                [coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)+params.coeff_desm(3,1)*coeffs(5)+params.coeff_desm(3,2)*coeffs(6)]]; %lvl 1=3
        case 'speed'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            regl = [[coeffs(1)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4),coeffs(2)+params.coeff_desm(3,1)*coeffs(5)+params.coeff_desm(3,2)*coeffs(6)]; ...
                [coeffs(1)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4),coeffs(2)+params.coeff_desm(1,1)*coeffs(5)+params.coeff_desm(1,2)*coeffs(6)]; ... %lvl 1
                [coeffs(1)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4),coeffs(2)+params.coeff_desm(2,1)*coeffs(5)+params.coeff_desm(2,2)*coeffs(6)]]; %lvl 1=3
        case 'speed:kin'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            % regl = [[coeffs(1),coeffs(2)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4)]; ...
            %     [coeffs(1),coeffs(2)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4)]; ... %lvl 1
            %     [coeffs(1),coeffs(2)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4)]]; %lvl 1=3
            regl = [coeffs(1),coeffs(2),coeffs(3),coeffs(4)];
        case 'speed:kin:grp'
            %-- (04/17/2025) JS, coeff organization for speed study OA vs YA
            % regl = [[coeffs(1),coeffs(2)+params.coeff_desm(3,1)*coeffs(3)+params.coeff_desm(3,2)*coeffs(4)]; ...
            %     [coeffs(1),coeffs(2)+params.coeff_desm(1,1)*coeffs(3)+params.coeff_desm(1,2)*coeffs(4)]; ... %lvl 1
            %     [coeffs(1),coeffs(2)+params.coeff_desm(2,1)*coeffs(3)+params.coeff_desm(2,2)*coeffs(4)]]; %lvl 1=3
            regl = coeffs;
    end

    %## ASSIGN STATS
    str = {'','',''};
    switch EXTRACT_STRUCT.str_option
        case 'simple'
            [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
        case 'information'
            [str] = get_inform_str(anv_p,tmp_stats.r2_c_int,regl);
    end
    % if EXTRACT_STRUCT.do_display_str
    %     % [str] = get_anv_str(anvs,tmp_stats.r2_c_int,fsqs);
    %     [str] = get_anv_str(anv_p,tmp_stats.r2_c_int,tmp_stats.f2_c_int);
    %     [str] = get_inform_str(anv_p,r2s,regl);
    % end

    %## RAW OUTPUT
    %--
    raw_values = struct('anova_p',anv_p, ...
        'anova_s',anv_stat, ...
        'coeffs',coeffs, ...
        'lme_ccs',{tmp_lme_cc}, ...
        'lme_pv',lme_pv, ...
        'lme_est',lme_est, ...
        'lme_stat',lme_stat, ...
        'fsq',fsqs, ...
        'rsq',[tmp_stats.r2_c_int,tmp_stats.r2_m_int]);

    %##  
    tssc = struct('var_type','continuous', ...
        'anova',anv_p(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',regl, ...
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',(0:5)*0.25,... % continuous predictors
        'order',{{}});
    tssg = struct('var_type','categorical', ...
        'anova',anv_p(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',[],... 
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',0,... % continuous predictors
        'order',{params.group_order});
    STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
        'cond_stats',tssc, ...
        'group_stats',tssg, ...
        'stats_char',struct('str',{str}, ...
            'offsets',params.str_offset, ...
            'font_size',params.str_font_size, ...
            'do_display',true));
end

end
%% ===================================================================== %%
function [str] = get_anv_str(anova_vals,r2s,efs)
    %##
    % if anova_vals(3) < 0.05 && anova_vals(3) > 0.01
    %     strg = '^{+}';
    % elseif anova_vals(3) <= 0.01 && anova_vals(3) > 0.001
    %     strg = '^{++}';
    % elseif anova_vals(3) <= 0.001
    %     strg = '^{+++}';
    % else
    %     strg = ''; %'^{ns}';
    % end
    % %--
    % if anova_vals(2) < 0.05 && anova_vals(2) > 0.01
    %     strs = '^{*}';
    % elseif anova_vals(2) <= 0.01 && anova_vals(2) > 0.001
    %     strs = '^{**}';
    % elseif anova_vals(2) <= 0.001
    %     strs = '^{***}';
    % else
    %     strs = ''; %'^{ns}';
    % end
    %##
    if anova_vals(3) < 0.05 && anova_vals(3) > 0.01
        strg = '+';
    elseif anova_vals(3) <= 0.01 && anova_vals(3) > 0.001
        strg = '++';
    elseif anova_vals(3) <= 0.001
        strg = '+++';
    else
        strg = ''; %'^{ns}';
    end
    %--
    if anova_vals(2) < 0.05 && anova_vals(2) > 0.01
        strs = '*';
    elseif anova_vals(2) <= 0.01 && anova_vals(2) > 0.001
        strs = '**';
    elseif anova_vals(2) <= 0.001
        strs = '***';
    else
        strs = ''; %'^{ns}';
    end
    %--
    if length(anova_vals) > 3
        % if anova_vals(4) < 0.05 && anova_vals(4) > 0.01
        %     stri = sprintf('^{%s}',char(8225));
        % elseif anova_vals(4) <= 0.01 && anova_vals(4) > 0.001
        %     stri = sprintf('^{%s%s}',char(8225),char(8225));
        % elseif anova_vals(4) <= 0.001
        %     stri = sprintf('^{%s%s%s}',char(8225),char(8225),char(8225));
        % else
        %     stri = ''; %'^{ns}';
        % end
        if anova_vals(4) < 0.05 && anova_vals(4) > 0.01
            stri = sprintf('%s',char(8225));
        elseif anova_vals(4) <= 0.01 && anova_vals(4) > 0.001
            stri = sprintf('%s%s',char(8225),char(8225));
        elseif anova_vals(4) <= 0.001
            stri = sprintf('%s%s%s',char(8225),char(8225),char(8225));
        else
            stri = ''; 
        end
        %-- make string
        % str = {sprintf('%sf_{s}^{2}=%1.2f   %sf_{g}^{2}=%1.2f   %sf_{s:g}^{2}=%1.2f\nR^2=%1.2f', ...
        %     strs,efs(1), ...
        %     strg,efs(2), ...
        %     stri,efs(3), ...
        %     r2s),'',''};
        %-- make string
        % str = {sprintf('%sR^{2}=%1.2f    f^{2}=%1.2f', ...
        %     strs, ...
        %     r2s, ...
        %     efs(1)),'',''};
        % str = {sprintf('%sR^{2}=%1.2f    f^{2}=%1.2f', ...
        %     stri, ...
        %     r2s, ...
        %     efs(1)),'',''};
        if ~isempty(strs)
            strs = sprintf('S^{%s}',strs);
        end
        if ~isempty(strg)
            strg = sprintf('G^{%s}',strg);
        end
        if ~isempty(stri)
            stri = sprintf('S:G^{%s}',stri);
        end
        sigs = {strs,strg,stri};
        sigs = sigs(~cellfun(@isempty,sigs));
        % str = {sprintf('%s\nR^{2}=%1.2f    f^{2}=%1.2f', ...
        %         strjoin(sigs,'  '), ...
        %         r2s, ...
        %         efs(1)),'',''};
        str = {sprintf('%s\nR^{2}=%1.2f', ...
                strjoin(sigs,'  '), ...
                r2s),'',''};
    else
        %-- make string
        % str = {sprintf('%sf_{s}^{2}=%1.2f   %sf_{g}^{2}=%1.2f\nR^2=%1.2f', ...
        %     strs,efs(1), ...
        %     strg,efs(2), ...
        %     r2s),'',''};
        %-- make string
        % str = {sprintf('%sR^{2}=%1.2f    f^{2}=%1.2f\n', ...
        %     strs, ...
        %     r2s, ...
        %     efs(1)),'',''};
        if ~isempty(strs)
            strs = sprintf('S^{%s}',strs);
        end
        if ~isempty(strg)
            strg = sprintf('G^{%s}',strg);
        end
        sigs = {strs,strg};
        sigs = sigs(~cellfun(@isempty,sigs));
        % str = {sprintf('%s\nR^{2}=%1.2f    f^{2}=%1.2f', ...
        %         strjoin(sigs,'  '), ...
        %         r2s, ...
        %         efs(1)),'',''};
        str = {sprintf('%s\nR^{2}=%1.2f', ...
                strjoin(sigs,'  '), ...
                r2s),'',''};
    end
    
    %-- make string
    % str = {sprintf('%sf_{s}^{2}=%1.2f   %sf_{g}^{2}=%1.2f   %sf_{s:g}^{2}=%1.2f\nR^2=%1.2f', ...
    %     strs,efs(1), ...
    %     strg,efs(2), ...
    %     stri,efs(3), ...
    %     r2s),'',''};
    %--
    % str = {sprintf('%sm_{ya}=%1.2f  m_{ohf}=%1.2f  m_{olf}=%1.2f\nR^2=%1.2f', ...
    %     stri,coeffs(2), ...
    %     coeffs(2)+coeffs(5), ...
    %     coeffs(2)+coeffs(6), ...
    %     tmp_stats.r2_c_int),'',''};
    % txt_sz = 9;
    % offs = [-0.11,-0.05];
    %--
    % str = {sprintf('%sy=(%1.1f)x+(%1.1f)\nR^2=%1.2f',stri,coeffs(2),coeffs(1),tmp_stats.r2_c_int), ...
    %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(5),coeffs(1)+coeffs(3)), ...
    %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(6),coeffs(1)+coeffs(4))};
end
function [str] = get_inform_str(anova_vals,r2s,coeffs)
    %##
    if anova_vals(3) < 0.05 && anova_vals(3) > 0.01
        strg = '+';
    elseif anova_vals(3) <= 0.01 && anova_vals(3) > 0.001
        strg = '++';
    elseif anova_vals(3) <= 0.001
        strg = '+++';
    else
        strg = ''; %'^{ns}';
    end
    %--
    if anova_vals(2) < 0.05 && anova_vals(2) > 0.01
        strs = '*';
    elseif anova_vals(2) <= 0.01 && anova_vals(2) > 0.001
        strs = '**';
    elseif anova_vals(2) <= 0.001
        strs = '***';
    else
        strs = ''; %'^{ns}';
    end
    %--
    coeffs = round(coeffs,2,'significant');
    fmt_num = '(%1.2f)';
    eq_str = {sprintf(['y_{YA}=',fmt_num,'x+',fmt_num],coeffs(1,2),coeffs(1,1)), ...
        sprintf(['y_{OHFA}=',fmt_num,'x+',fmt_num],coeffs(2,2),coeffs(2,1)), ...
        sprintf(['y_{OLFA}=',fmt_num,'x+',fmt_num],coeffs(3,2),coeffs(3,1))};
    %--
    if length(anova_vals) > 3
        if anova_vals(4) < 0.05 && anova_vals(4) > 0.01
            stri = sprintf('%s',char(8225));
        elseif anova_vals(4) <= 0.01 && anova_vals(4) > 0.001
            stri = sprintf('%s%s',char(8225),char(8225));
        elseif anova_vals(4) <= 0.001
            stri = sprintf('%s%s%s',char(8225),char(8225),char(8225));
        else
            stri = ''; 
        end
        if ~isempty(strs)
            strs = sprintf('S^{%s}',strs);
        end
        if ~isempty(strg)
            strg = sprintf('G^{%s}',strg);
        end
        if ~isempty(stri)
            stri = sprintf('S:G^{%s}',stri);
        end
        sigs = {strs,strg,stri};
        sigs = sigs(~cellfun(@isempty,sigs));
        %--
        % str = {sprintf('%s\nR^{2}=%1.2f\n%s', ...
        %         strjoin(sigs,'  '), ...
        %         r2s,eq_str{1}),eq_str{2},eq_str{3}};
        str = {sprintf('%s\nR^{2}=%1.2f\n%s', ...
                strjoin(sigs,'  '), ...
                r2s,strjoin(eq_str,'  ')),'',''};
    else
        if ~isempty(strs)
            strs = sprintf('S^{%s}',strs);
        end
        if ~isempty(strg)
            strg = sprintf('G^{%s}',strg);
        end
        sigs = {strs,strg};
        sigs = sigs(~cellfun(@isempty,sigs));
        % str = {sprintf('%s\nR^{2}=%1.2f\n%s', ...
        %         strjoin(sigs,'  '), ...
        %         r2s,eq_str{1}),eq_str{2},eq_str{3}};
        str = {sprintf('%s\nR^{2}=%1.2f\n%s', ...
                strjoin(sigs,'  '), ...
                r2s,strjoin(eq_str,'  ')),'',''};
    end
    
end

function tmp_var = tbl2num(match_chars,tbl_chars,tbl,col_char)
    %-- 
    tmp = tbl.(col_char);
    dblv = cellfun(@(x) double(string(x)),strsplit(tmp{1},','));
    tmp_var = zeros(length(match_chars),1);
    for cc = 1:length(match_chars)
        ind = strcmp(match_chars{cc},tbl_chars);
        if ~isempty(ind)
            tmp_var(cc) = dblv(ind);              
        else
            fprintf("Coefficient %s not found.\n",match_chars{cc})
        end         
    end
end
%% =======================================================================
%## POTENTIAL CODE FOR LATER DEV
%-- model coefficients
% coeff_chars = strcmp(params.group_chars,'(Intercept)');
% coeff_chars = params.group_chars(~coeff_chars);    
% % params.coeff_chars_unmix = {'','group_char1','group_char2';'OHFA','OLFA','YA'};
% params.coeff_chars_unmix = {'(Intercept)','speed_cond_num','group_char1','group_char2','speed_cond_num:group_char1','speed_cond_num:group_char2'; ...
%     'OHFA','OHFA','OLFA','YA','OLFA','YA'};
% % coeffs = zeros(length(params.coeff_chars_int),1);
% coeffs = zeros(length(params.group_chars),length(params.coeff_chars_int));
% for cc = 1:length(params.coeff_chars_int)
%     ind = strcmp(params.coeff_chars_int{cc},tmp_cc);
%     %-- unmix
%     % unind = cellfun(@(x) contains(tmp_cc{ind},x) && ~isempty(x),params.coeff_chars_unmix(1,:));        
%     % bind = cellfun(@isempty,params.coeff_chars_unmix(1,:));
%     %--
%     unind = cellfun(@(x) strcmp(tmp_cc{ind},x),params.coeff_chars_unmix(1,:));
%     if any(unind) && ~isempty(ind)
%         unindun = strcmp(params.coeff_chars_unmix{2,unind},coeff_chars);
%         % tmp = tmp_coeffs(ind);
%         coeffs(unindun,cc) = tmp_coeffs(ind);
%     % elseif ~isempty(ind)
%     %     unindun = strcmp(params.coeff_chars_unmix{2,bind},coeff_chars);
%     %     coeffs(1,cc) = tmp_coeffs(ind); 
%     else
%         error("Coefficient %s not found.\n",params.coeff_chars_int{cc});            
%     end        
% end
