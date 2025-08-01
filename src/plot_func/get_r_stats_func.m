function [ranef,raw_values,tmp_stats] = get_r_stats_func(rstats_table,cl_n,varargin)
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
%--
DEF_EXTRACT_STRUCT = struct( ...
    'model_char', 'speed_group_intact_all', ...
    'group_char', 'all', ...
    'anv_chars', {ANV_CHARS_INT}, ...
    'coeff_chars', {COEFF_CHARS_INT}, ...
    'eeg_measure', {''}, ...
    'kin_measure', {'none'}, ...
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

%% ===================================================================== %%
%## EXTRACT STATS INFO
if ~isfield(params,'kin_measure')
    params.kin_measure = [];
end
if isempty(params.kin_measure)
    % tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
    %     strcmp(rstats_table.model_char,params.model_char) &...
    %     strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
    %     strcmp(rstats_table.group_char,params.group_char);
    tmp_stats = rstats_table.cluster_n==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
else
    % tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
    %     strcmp(rstats_table.model_char,params.model_char) &...
    %     strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
    %     strcmp(rstats_table.kinematic_char,params.kin_measure) &...
    %     strcmp(rstats_table.group_char,params.group_char);
    tmp_stats = rstats_table.cluster_n==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.kinematic_char,params.kin_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
end
tmp_stats = rstats_table(tmp_stats,:);

%## EXTRACT CHARACTERS
tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');

%## EXTRACT NUMBERS
%-- rndm inter
tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
ranef = struct('char',tmp_ranef_char, ...
    'int',tmp_ranef_n);    

%-- lme params
try
    tmp_lme_cc = strsplit(tmp_stats.lme_cc{1},',');
    lme_est = tbl2num(params.coeff_chars,tmp_lme_cc,tmp_stats,'lme_est');
    lme_stat = tbl2num(params.coeff_chars,tmp_lme_cc,tmp_stats,'lme_stat');
    lme_pv = tbl2num(params.coeff_chars,tmp_lme_cc,tmp_stats,'lme_pval');
catch
    fprintf('lme_cc doesn''t exist.\n');
    tmp_lme_cc = {''};
    lme_est = [];
    lme_stat = [];
    lme_pv = [];
end
coeffs = tbl2num(params.coeff_chars,tmp_cc,tmp_stats,'coeffs');
%-- anova params
anv_p = tbl2num(params.anv_chars,tmp_ac,tmp_stats,'anv_pvals');
anv_stat = tbl2num(params.anv_chars,tmp_ac,tmp_stats,'anv_stats');
%-- emm params
if ~isempty(tmp_fsq_chars)
    fsqs = tbl2num(tmp_fsq_chars,tmp_fsq_chars,tmp_stats,'fsq_vals');
else
    fsqs = [];
end
if ~isempty(tmp_ci_chars)
    tmp_ci_lwr = tbl2num(tmp_ci_chars,tmp_ci_chars,tmp_stats,'confint_lwr');
    tmp_ci_upr = tbl2num(tmp_ci_chars,tmp_ci_chars,tmp_stats,'confint_upr');
    cis = zeros(length(tmp_ci_chars),1,2);
    for cc = 1:length(tmp_ci_chars)
        cis(cc,1,:) = [tmp_ci_lwr(cc),tmp_ci_upr(cc)];
    end
else
    cis = [];
end
%-- confidence intervals
% if ~all(cellfun(@isempty,tmp_ci_chars))
%     chkd = false;
%     CONFINT_STRUCT = struct('do_display',chkd, ...
%             'y_bnds',cis, ...
%             'x_vals',repmat(params.ci_bar_xpos,[3,1]), ...
%             'errbar_struct',struct('line_specs',{params.ci_bar_linespecs}, ...
%                 'err_bar_width',params.ci_bar_width));
% end
%-- raw values
raw_values = struct('anova_ccs',{params.anv_chars}, ... %{tmp_ac}, ...
    'anova_pv',anv_p, ...
    'anova_stat',anv_stat, ...
    'coeffs',coeffs, ...
    'lme_ccs',{params.coeff_chars}, ... %{tmp_lme_cc}, ...
    'lme_pv',lme_pv, ...
    'lme_est',lme_est, ...
    'lme_stat',lme_stat, ...
    'cis_chars',{tmp_ci_chars}, ...
    'cis_lwrupr',cis, ...
    'fsq_chars',{tmp_fsq_chars}, ...
    'fsq',fsqs, ...
    'rsq',[tmp_stats.r2_c_int,tmp_stats.r2_m_int]);

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
