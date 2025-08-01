function [ax] = group_violin(table_in,measure_char,cond_char,group_char,varargin)
%GROUP_VIOLIN Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%- FONTS
% ax = [];
SCATTER_MARK_SIZE = 5;
STD_CUTOFF = 3;
DEF_VIOLIN_STRUCT = struct('Width',0.05,...
            'ShowWhiskers',false,...
            'ShowNotches',false,...
            'ShowBox',true,...
            'ShowMedian',true,...
            'Bandwidth',0.1,...
            'QuartileStyle','shadow',...
            'HalfViolin','full',...
            'DataStyle','scatter',...
            'MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],...
            'ViolinAlpha',{{0.2 0.3}},...
            'vio_edge_alpha',0.9, ...
            'vio_edge_lw',0.5, ...
            'do_plot_outlier_marks',true,...
            'use_raw_bandwidth',false);
violin_params_keep = {'Width','ShowWhiskers','ShowNotches','ShowBox', ...
    'ShowMedian','Bandwidth','QuartileStyle','HalfViolin', ...
    'DataStyle','MarkerSize','EdgeColor','ViolinAlpha'};
label_struct = struct('Units','normalized',...
    'FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
group_anno_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center',...
    'Units','data');
ax_struct = struct('box','off', ...
    'LineWidth',1, ...
    'FontWeight','normal', ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'OuterPosition',[0 0 1 1], ...
    'Position',[0.3,0.3,0.5,0.5]);
title_props_struct = struct('Units','normalized',...
    'FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
stats_char_struct = struct('FontSize',9,...
    'FontName','Arial',...
    'FontWeight','bold', ...
    'Units','normalized');
reg_line_struct = struct('LineStyle','-', ...
    'Color','k', ...
    'LineWidth',2);
DEF_PLOT_STRUCT = struct('color_map',linspecer(length(unique(table_in.(cond_char)))),...
    'cond_labels',{cellstr(string(unique(table_in.(cond_char))))},...
    'cond_offsets',linspace(-0.3,0.3,length(unique(table_in.(cond_char)))),...
    'do_group_labels',true, ... % group props
    'do_combine_groups',false,...
    'group_labels',{cellstr(string(unique(table_in.(group_char))))},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_label_offset',[0,-0.285],...
    'group_label_props',group_anno_struct, ...
    'y_label',{measure_char},... % y-label props
    'ytick',[], ... 
    'ytick_labs',{{''}}, ...
    'y_label_props',label_struct, ...
    'y_label_offset',[0,0], ...
    'ylim',[],...
    'x_label',{cond_char},... % x-label props
    'x_label_props',label_struct, ...
    'x_label_offset',[0,-0.1],...
    'xlim',[],...
    'xtick',[], ...
    'xtick_labs',{{''}}, ...
    'xtick_angle',75,...
    'title',{{''}},... % title props
    'title_props',title_props_struct, ...
    'ax_props',ax_struct, ... % ax props
    'stats_char_props',stats_char_struct, ... % stats & scatter props
    'reg_line_props',reg_line_struct, ...
    'scatter_props',struct(...
        'SizeData',4, ...
        'Marker','*', ...
        'CData',[0,0,0], ...
        'jitter','on', ...
        'jitterAmount', 0.05));
%## STATISTICS STRUCTURE    
tss = struct('var_type','numeric', ...
    'anova',1, ...
    'multc_pvals',[],...
    'multc_pvals_pairs',[],...
    'regress_line',[],... 
    'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
    'regress_xvals',0,... % continous indepedent variable
    'order',{{}});
cond_stats = tss;
group_stats = tss;
% cond_stats.order = cellstr(string(unique(table_in.(cond_char))));
% group_stats.order = cellstr(string(unique(table_in.(group_char))));
DEF_STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
    'cond_stats',cond_stats, ...
    'group_stats',group_stats, ...
    'stats_char',struct('str',{{}}, ...
        'offsets',[-0.1,0], ...        
        'font_size',9, ...
        'do_display',true));
%## STATS VISUALIZATION STRUCTURES
DEF_BRACKET_STRUCT = struct('sig_sign','+',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
DEF_SIGLINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_offset_x',0,...
    'sig_offset_y',0);
DEF_CONFINT_STRUCT = struct('do_display',false, ...
    'y_bnds',[], ... % size: [group, condition, 2]
    'x_vals',[], ... % size: [group, condition]
    'errbar_struct',struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
        'err_bar_width',1));
%- STATS_STRUCT EXAMPLE:
%## CHECK FUNCTIONS
isnumscalar = @(x) (isnumeric(x)&isscalar(x));
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'table_in',@istable);
addRequired(p,'measure_char',@ischar);
addRequired(p,'cond_char',@ischar);
addRequired(p,'group_char',@ischar);
%## OPTIONAL
addOptional(p,'parent_axes',[],@(x) isobject(x) || isempty(x))
%## PARAMETER
addParameter(p,'VIOLIN_STRUCT',DEF_VIOLIN_STRUCT,@(x) validate_struct(x,DEF_VIOLIN_STRUCT));
addParameter(p,'PLOT_STRUCT',DEF_PLOT_STRUCT,@(x) validate_struct(x,DEF_PLOT_STRUCT));
addParameter(p,'STATS_STRUCT',DEF_STATS_STRUCT,@(x) validate_struct(x,DEF_STATS_STRUCT));
addParameter(p,'BRACKET_STRUCT',DEF_BRACKET_STRUCT,@(x) validate_struct(x,DEF_BRACKET_STRUCT));
addParameter(p,'SIGLINE_STRUCT',DEF_SIGLINE_STRUCT,@(x) validate_struct(x,DEF_SIGLINE_STRUCT))
addParameter(p,'CONFINT_STRUCT',DEF_CONFINT_STRUCT,@(x) validate_struct(x,DEF_CONFINT_STRUCT))
%##
parse(p,table_in,measure_char,cond_char,group_char,varargin{:});
%## SET DEFAULTS
ax = p.Results.parent_axes;
%-- extract vars
VIOLIN_STRUCT = p.Results.VIOLIN_STRUCT;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
STATS_STRUCT = p.Results.STATS_STRUCT;
BRACKET_STRUCT = p.Results.BRACKET_STRUCT;
SIGLINE_STRUCT = p.Results.SIGLINE_STRUCT;
CONFINT_STRUCT = p.Results.CONFINT_STRUCT;
%-- set defaults
VIOLIN_STRUCT = set_defaults_struct(VIOLIN_STRUCT,DEF_VIOLIN_STRUCT);
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEF_PLOT_STRUCT);
STATS_STRUCT = set_defaults_struct(STATS_STRUCT,DEF_STATS_STRUCT);
BRACKET_STRUCT = set_defaults_struct(BRACKET_STRUCT,DEF_BRACKET_STRUCT);
SIGLINE_STRUCT = set_defaults_struct(SIGLINE_STRUCT,DEF_SIGLINE_STRUCT); 
CONFINT_STRUCT = set_defaults_struct(CONFINT_STRUCT,DEF_CONFINT_STRUCT); 
SIGLINE_STRUCT.sig_levels = STATS_STRUCT.sig_levels;

%## CONVERT STRUCT CELL ARGS
x_label_props = struct2args(PLOT_STRUCT.x_label_props);
y_label_props = struct2args(PLOT_STRUCT.y_label_props);
group_label_props = struct2args(PLOT_STRUCT.group_label_props);
title_props = struct2args(PLOT_STRUCT.title_props);
ax_props = struct2args(PLOT_STRUCT.ax_props);
stats_char_props = struct2args(PLOT_STRUCT.stats_char_props);
reg_line_props = struct2args(PLOT_STRUCT.reg_line_props);
scatter_props = struct2args(PLOT_STRUCT.scatter_props);
%% ===================================================================== %%
%##
group_stats = STATS_STRUCT.group_stats;
cond_stats = STATS_STRUCT.cond_stats;
stats_char = STATS_STRUCT.stats_char;
table_tmp = table_in;
%-- type checks
ccchk = false;
cnchk = false;
gcchk = false;
gnchk = false;
%-- cond chk
switch cond_stats.var_type
    case {'numeric','num','continuous','cont'}
        cnchk = true;
    case {'categorical','cat','factor','fac'}
        ccchk = true;
    otherwise
        error('Enter a valid STATS_STRUCT.cond_stats.var_type: ''numeric'',''continuous''');
end
%-- group chk
switch group_stats.var_type
    case {'numeric','num','continuous','cont'}
        gnchk = true;
    case {'categorical','cat','factor','fac'}
        gcchk = true;
    otherwise
        error('Enter a valid STATS_STRUCT.cond_stats.var_type: ''numeric'',''continuous''');
end

%## SET GROUPS
%-- set group order
if any(isempty(group_stats.order)) || any(isundefined(group_stats.order))
    groups = unique(table_tmp.(group_char));
else
    warning('STATS_STRUCT.group_stats.order is being used, make sure to change PLOT_STRUCT.group_labels to match...')
    groups = group_stats.order;
    if ~iscategorical(groups)
        groups = categorical(groups);
    end
end
%-- combine groups if so
if PLOT_STRUCT.do_combine_groups
    table_tmp.(group_char) = categorical(ones(height(table_tmp),1));
    groups = unique(table_tmp.(group_char));
else
    fprintf('Plotting groups seperately...\n')
end
%-- set conditions
if any(isempty(cond_stats.order)) || any(isundefined(cond_stats.order))
    conds = unique(table_tmp.(cond_char));
else
    warning('STATS_STRUCT.cond_stats.order is being used, make sure to change PLOT_STRUCT.cond_labels to match...')
    conds = cond_stats.order;
    if ~iscategorical(conds)
        conds = categorical(conds);
    end
end
%% MAIN VIOLIN PLOTTING LOO ============================================ %%
hold on;
cnt = 1;
g_o = 1; 
xticks = [];
xtick_labs = {};
gticks = [];
y_stats = zeros(1,30);
y_cnt = 1;
%-- group primary, cond secondary
ymaxs = zeros(length(groups),length(conds));
sigline_ymax = zeros(length(groups),length(conds));
violins = cell(length(groups)*length(conds),1);
if ~isempty(PLOT_STRUCT.group_offsets)
    g_offset = PLOT_STRUCT.group_offsets(g_o);
    g_o = g_o + 1;
else
    g_offset = 0;
end
tmp_bandwidth = VIOLIN_STRUCT.Bandwidth;
stats_design_mat = zeros(length(groups),length(conds),2);
for i=1:length(groups)
    for j=1:length(conds)
        g_i = groups(i);
        c_i = conds(j);
        %-- stats design mat
        stats_design_mat(i,j,1) = i;
        stats_design_mat(i,j,2) = j;
        %-- typing checks
        if iscell(table_tmp.(cond_char))
            table_tmp.(cond_char) = categorical(table_tmp.(cond_char));
        end
        if iscell(table_tmp.(group_char))
            table_tmp.(group_char) = categorical(table_tmp.(group_char));
        end
        inds = table_tmp.(cond_char)==c_i & table_tmp.(group_char)==g_i;
        %-- outliers
        tt = table_tmp(inds,:);
        vals = tt.(measure_char);
        sigline_ymax(i,j) = max(tt.(measure_char));
        chk = isnan(vals);
        %-- remove NaNs from vals
        if any(chk)
            % if ~isempty(STATS_STRUCT.subject_char)
            %     fprintf('\nRemoving NaN''s from data for subjects: \n');
            %     fprintf('%s, ',tt.(STATS_STRUCT.subject_char)(logical(chk)))
            % end
            vals = vals(~chk);
        end
        %-- calculate mean and std and crop out data beyond limits
        if VIOLIN_STRUCT.do_plot_outlier_marks
            data_mean = mean(vals);
            data_std = std(vals);
            ind_crop = tt.(measure_char)>data_mean-STD_CUTOFF*data_std & tt.(measure_char)<data_mean+STD_CUTOFF*data_std;
            outliers = tt(~ind_crop,:);
            tt = tt(ind_crop,:);
        else
            outliers = [];
        end
        %-- final table in
        offset = PLOT_STRUCT.cond_offsets(j)+g_offset;
        %-- check input types
        if ~isnumscalar(g_i)
            tmp_gi = i; 
        else
            tmp_gi = g_i;
        end
        %## MAIN MEASURES
        if ~isempty(tt)
            ymaxs(i,j) = max(tt.(measure_char));

            %## VIOLIN STRUCT EDITS
            %-- plot main violin for each condition & group
            if ~VIOLIN_STRUCT.use_raw_bandwidth
                VIOLIN_STRUCT.Bandwidth = range(vals)*tmp_bandwidth;
            end
            %-- clean up violin_params
            VIOLIN_PARAMS = struct2args(VIOLIN_STRUCT);
            fns = fieldnames(VIOLIN_STRUCT);
            inds = find(cellfun(@(x) any(strcmp(x,violin_params_keep)),fns));
            keep_inds = zeros(length(VIOLIN_PARAMS),1);
            for ff = 1:length(inds)
                ind = find(strcmp(fns(inds(ff)),VIOLIN_PARAMS));
                keep_inds(ind) = 1;
                keep_inds(ind+1) = 1;
            end
            VIOLIN_PARAMS = VIOLIN_PARAMS(logical(keep_inds));

            %## MAKE VIOLIN PLOT
            % violins(cnt) = Violin({table_in.(measure_char)},tmp_gi,VIOLIN_PARAMS{:});
            violins{cnt} = Violin({tt.(measure_char)},tmp_gi,VIOLIN_PARAMS{:});
            violins{cnt}.ScatterPlot.XData      = violins{cnt}.ScatterPlot.XData+offset;
            violins{cnt}.ViolinPlot.XData       = violins{cnt}.ViolinPlot.XData+offset;
            violins{cnt}.ViolinPlot.EdgeAlpha   = VIOLIN_STRUCT.vio_edge_alpha;
            violins{cnt}.ViolinPlot.LineWidth   = VIOLIN_STRUCT.vio_edge_lw; 
            try
                violins{cnt}.WhiskerPlot.XData      = violins{cnt}.WhiskerPlot.XData+offset;
                
            catch e
                fprintf('Error. %s\n%s\n',e.identifier,getReport(e));
            end
            violins{cnt}.MedianPlot.XData       = violins{cnt}.MedianPlot.XData+offset;
            violins{cnt}.NotchPlots(1).XData    = violins{cnt}.NotchPlots(1).XData+offset;
            violins{cnt}.NotchPlots(2).XData    = violins{cnt}.NotchPlots(2).XData+offset;
            violins{cnt}.MeanPlot.XData         = violins{cnt}.MeanPlot.XData+offset;            
            violins{cnt}.ViolinPlotQ.XData      = violins{cnt}.ViolinPlotQ.XData+offset;
            %-- set some props
            violins{cnt}.MeanPlot.LineWidth     = VIOLIN_STRUCT.vio_edge_lw;
            violins{cnt}.ViolinPlotQ.LineWidth  = VIOLIN_STRUCT.vio_edge_lw;
            % violins{cnt}.ViolinPlotQ.LineStyle  = '-';
            violins{cnt}.WhiskerPlot.LineWidth  = VIOLIN_STRUCT.vio_edge_lw; 
            violins{cnt}.ViolinColor            = {PLOT_STRUCT.color_map(j,:)};
            %-- get xticks and labels
            xticks = [xticks,tmp_gi+offset];
            if iscell(PLOT_STRUCT.cond_labels)
                xtick_labs = [xtick_labs,{sprintf('%s',PLOT_STRUCT.cond_labels{j})}];
            else
                xtick_labs = [xtick_labs,{string(PLOT_STRUCT.cond_labels(j))}];
            end
            %## OUTLIERS
            if ~isempty(outliers)
                g_ind = (outliers.(cond_char)==c_i) & (outliers.(group_char)==g_i);
                vals_y = outliers.(measure_char)(g_ind);
                if ~isempty(vals_y)
                    in_x = repmat(xticks(cnt),size(vals_y));
                    ss = scatter(in_x,vals_y,12,'k*');
                    set(ss,scatter_props{:});
                    y_stats(y_cnt) = gety(ss);
                    y_cnt = y_cnt + 1;
                end
            else
                fprintf('Condition %s & Group %s does not have outliers\n',string(c_i),string(g_i))
            end
        else
            fprintf('Condition %s & Group %s does not have entries\n',string(c_i),string(g_i))
        end        
        cnt = cnt+1;
    end
    %-- handle groupings offsets
    gticks(i) = mean(xticks((cnt-length(conds)):(cnt-1)));
    if ~isempty(PLOT_STRUCT.group_offsets) & i < length(groups)
        g_offset = PLOT_STRUCT.group_offsets(g_o);
        g_o = g_o + 1;
    else
        g_offset = g_offset + 0.2;
    end
end
hold on;

%## DISPLAY STATISTICS
hold_xlim = get(gca,'xlim');
% set_y = true;
stats_displayed = false;
tmp_bracket_struct = BRACKET_STRUCT;
tmp_sigline_struct = SIGLINE_STRUCT;
tmp_bracket_struct.sig_levels = STATS_STRUCT.sig_levels;
tmp_sigline_struct.sig_levels = STATS_STRUCT.sig_levels;
%## CATEGORICAL GROUP STATISTICSv
% cnt_g = 1;
%-- check that stats exist
chk = ~isempty(group_stats.anova) && ...
    group_stats.anova < STATS_STRUCT.sig_levels(1) && ...
    gcchk;
if chk
    % for i=1:length(groups)
    for i = 1:length(group_stats.multc_pvals_pairs)
        % if STATS_STRUCT.pvals_grp{i} < STATS_STRUCT.sig_levels(1)
        if group_stats.multc_pvals(i) < STATS_STRUCT.sig_levels(1)
            % g1 = STATS_STRUCT.pvals_grp_pairs{i}(1,1);
            % g2 = STATS_STRUCT.pvals_grp_pairs{i}(1,2);
            v1 = group_stats.multc_pvals_pairs(i,1);
            v2 = group_stats.multc_pvals_pairs(i,2);
            if v1 ~= 1
                bx1 = [(v1-1)*length(conds)+1,(v1-1)*length(conds)+length(conds)];
            else
                bx1 = [v1,v1+length(conds)-1];
            end
            if v2 ~= 1
                bx2 = [(v2-1)*length(conds)+1,(v2-1)*length(conds)+length(conds)];
            else
                bx2 = [v2,v2+length(conds)-1];
            end
            tmpy = [];
            for tt = linspace(bx1(1),bx1(2),bx1(2)-bx1(1)+1)
                tmpy = [tmpy, violins{tt}.ScatterPlot.YData];
            end
            by1 = [max(tmpy)*1.025,max(tmpy)*1.05];
            tmpy = [];
            for tt = linspace(bx2(1),bx2(2),bx2(2)-bx2(1)+1)
                tmpy = [tmpy, violins{tt}.ScatterPlot.YData];
            end
            by2 = [max(tmpy)*1.025,max(tmpy)*1.05];
            %-
            for tt = 1:length(bx1)
                bx1(tt) = violins{bx1(tt)}.MedianPlot.XData;
                bx2(tt) = violins{bx2(tt)}.MedianPlot.XData;
            end
            %- bracket 1
            bracket_1 = [[bx1(1),bx1(1),bx1(2),bx1(2)];...
                [by1(1),by1(2),by1(2),by1(1)]];
            %- bracket 2
            bracket_2 = [[bx2(1),bx2(1),bx2(2),bx2(2)];...
                [by2(1),by2(2),by2(2),by2(1)]];
            if length(BRACKET_STRUCT.bracket_offset_y_lower) > 1
                tmp_bracket_struct.bracket_offset_y_lower = BRACKET_STRUCT.bracket_offset_y_lower(i);
            end
            if length(BRACKET_STRUCT.conn_offset_y_upper) > 1
                tmp_bracket_struct.conn_offset_y_upper = BRACKET_STRUCT.conn_offset_y_upper(i);
            end
            if length(BRACKET_STRUCT.bracket_offset_y_upper) > 1
                tmp_bracket_struct.bracket_offset_y_upper = BRACKET_STRUCT.sig_offset_y(i);
            end
            if length(BRACKET_STRUCT.sig_offset_x) > 1
                tmp_bracket_struct.sig_offset_x = BRACKET_STRUCT.sig_offset_y(i);
            end
            if length(BRACKET_STRUCT.sig_offset_y) > 1
                tmp_bracket_struct.sig_offset_y = BRACKET_STRUCT.sig_offset_y(i);
            end
            pp = cus_sigbracket(gca,group_stats.multc_pvals(i),bracket_1,bracket_2,...
                'LINE_STRUCT',tmp_bracket_struct);
            y_stats(y_cnt) = gety(pp(3));
            y_cnt = y_cnt + 1;
        end
    end
end

%## CATEGORICAL CONDITION STATISTICS
%-- check that stats exist
chk = ~isempty(cond_stats.anova) && ...
    cond_stats.anova < STATS_STRUCT.sig_levels(1) && ...
    ccchk;
% sig_gap_y = [];
sig_gap_y = SIGLINE_STRUCT.sig_offset_y;
if chk
    for mc_i = 1:size(cond_stats.multc_pvals_pairs,1)
        %-- array with each comparison as a 2 set numbers of [group,condition]
        % where the 2nd dimension of the array is the comparion idx and the 3rd
        % dimension is the [group,condition] definition, and 1st dimension
        % is each comparison performed.
        ppair = squeeze(cond_stats.multc_pvals_pairs(mc_i,:,:));
        %-- array corresponding with multc_pvals_pairs at index pi
        pval = cond_stats.multc_pvals(mc_i);
        %-- get violin indices
        gsdm = stats_design_mat(:,:,1);
        csdm = stats_design_mat(:,:,2);
        v1 = find(gsdm == ppair(1,1) & csdm == ppair(1,2));
        v2 = find(gsdm == ppair(2,1) & csdm == ppair(2,2));
        %--
        % fpritnf('Plotting multiple comparisons for %s:%s',)
        xx(1) = violins{v1}.MedianPlot.XData; %#ok<FNDSB>
        xx(2) = violins{v2}.MedianPlot.XData; %#ok<FNDSB>
        if length(SIGLINE_STRUCT.sig_offset_x) > 1
            tmp_sigline_struct.sig_offset_x = SIGLINE_STRUCT.sig_offset_x(mc_i);
        end
        % if length(SIGLINE_STRUCT.sig_offset_y) > 1 && ~isempty(sig_gap_y)
        %     tmp_sigline_struct.sig_offset_y = SIGLINE_STRUCT.sig_offset_y(mc_i);
        % else 
        %     tmp_sigline_struct.sig_offset_y  = sig_gap_y;
        % end
        tmp_sigline_struct.sig_offset_y  = sig_gap_y;
        if length(SIGLINE_STRUCT.conn_y) > 1
            tmp_sigline_struct.conn_y = SIGLINE_STRUCT.conn_y(mc_i,:);
        end
        if length(SIGLINE_STRUCT.conn_offset_y) > 1
            tmp_sigline_struct.conn_offset_y = SIGLINE_STRUCT.conn_offset_y(mc_i);
        end
        [pp,sig_gap_y] = cus_sigline(gca,pval,xx,...
            'LINE_STRUCT',tmp_sigline_struct);
        y_stats(y_cnt) = gety(pp);
        y_cnt = y_cnt + 1;
    end
end

%## CONDITION REGRESSION-CONTINUOUS STATISTICS 
cnt_g = 1;
%-- check that stats exist
chk = cond_stats.anova < STATS_STRUCT.sig_levels(1) && ...
    ~isempty(cond_stats.regress_line) && ...
    cnchk;
if chk
    for i = 1:size(CONFINT_STRUCT.x_vals,1)
        %- plot line of best fit        
        if strcmp(cond_stats.line_type,'best_fit')
            m = cond_stats.regress_line(i,2);
            b = cond_stats.regress_line(i,1);
            %--
            xx = xticks(cnt_g:cnt_g+length(conds)-1);
            incr = xx(2)-xx(1);
            xx1 = [xx(1)-incr, xx];
            xx2 = [xx1, xx1(end)+incr];
            yy = cond_stats.regress_xvals*m + b;
            plot(xx2,yy,reg_line_props{:}); 
        %-
        elseif strcmp(cond_stats.line_type,'means')
            xx = xticks(cnt_g:cnt_g+length(conds)-1);
            yy = cond_stats.regress_line(i,:);
            plot(xx,yy,reg_line_props{:}) 
        end
        cnt_g = cnt_g + length(conds);
    end
end

%## CONFIDENCE INTERVALS
%-- check that stats exist
chk = CONFINT_STRUCT.do_display && ...
    ((cond_stats.anova < STATS_STRUCT.sig_levels(1) && ...
    ccchk) || (group_stats.anova < STATS_STRUCT.sig_levels(1) && ...
    gcchk));
%-- loop vars
cnt_g = 1;
tmp_errbar_struct = CONFINT_STRUCT.errbar_struct;
if chk
    for i = 1:size(CONFINT_STRUCT.x_vals,1)
        for j = 1:size(CONFINT_STRUCT.x_vals,2)
            %-- grab error bar positions
            tmp = xticks(cnt_g:cnt_g+length(conds)-1);                        
            xl = tmp(end)-tmp(1);
            fxl = CONFINT_STRUCT.errbar_struct.err_bar_width;
            if xl == 0
                xl = get(ax,'XLim');
                xl = (xl(2) - xl(1))/(size(CONFINT_STRUCT.x_vals,1)*size(CONFINT_STRUCT.x_vals,2));
            end
            tmp_errbar_struct.err_bar_width = fxl*xl;
            %-- shift relative positions
            if size(CONFINT_STRUCT.x_vals,2) > 1
                tt = tmp(j);
            else
                tt = mean(tmp);
            end
            x_vals = tt+CONFINT_STRUCT.x_vals(i,j);
            y_bnds = squeeze(CONFINT_STRUCT.y_bnds(i,j,:));
            %-- plot
            [~] = cus_err_bar(gca,x_vals,y_bnds, ...
                'ERRBAR_STRUCT',tmp_errbar_struct);
        end
        cnt_g = cnt_g + length(conds);
    end        
end

%## STATISTICS LABEL DISPLAY
%-- check that stats exist
chk = cond_stats.anova < STATS_STRUCT.sig_levels(1) || ...
    group_stats.anova < STATS_STRUCT.sig_levels(1);
if chk
    for i = 1:length(groups)
    %## STATS CHARACTER
        x_txt = PLOT_STRUCT.group_offsets(i)+stats_char.offsets(1);
        y_txt = 0.9+stats_char.offsets(2);
        if stats_char.do_display && ~stats_displayed
            if length(stats_char.str) ~= length(groups)
                error('ERROR. stats_char must be length(unique(table.group_char)).')
            end
            str = stats_char.str{i};
            tx = text(x_txt,y_txt,str,...
                stats_char_props{:});
            y_stats(y_cnt) = gety(tx);
            y_cnt = y_cnt + 1;
        end
    end
end

%## SET AXES PROPERTIES (AESTHETICS, LABELS, TICKS)
if isempty(ax)
    ax = gca;
end
set(ax,ax_props{:});

%## SET XLABEL
xh = xlabel(ax,PLOT_STRUCT.x_label, ...
    x_label_props{:});
%-- set pos
tp = get(xh,'Position');
tp(1,2) = tp(1,2)+PLOT_STRUCT.x_label_offset(2);
tp(1,1) = tp(1,1)+PLOT_STRUCT.x_label_offset(1);
set(xh,'Position',tp);
%-- xticks
if ~isempty(PLOT_STRUCT.xlim)
    set(ax,'XLim',PLOT_STRUCT.xlim)
else
    set(ax,'XLim',hold_xlim)
end
if ~isempty(PLOT_STRUCT.xtick)
    set(ax,'XTick',sort(PLOT_STRUCT.xtick));
else
    set(ax,'XTick',sort(xticks));
end
if ~all(cellfun(@isempty,PLOT_STRUCT.xtick_labs))
    set(ax,'XTickLabel',PLOT_STRUCT.xtick_labs);
else
    set(ax,'XTickLabel',xtick_labs);
end
xtickangle(PLOT_STRUCT.xtick_angle);

%## SET YLABEL
yh = ylabel(ax,PLOT_STRUCT.y_label, ...
    y_label_props{:});
%-- set pos
tp = get(yh,'Position');
tp(1,2) = tp(1,2)+PLOT_STRUCT.y_label_offset(2);
tp(1,1) = tp(1,1)+PLOT_STRUCT.y_label_offset(1);
set(yh,'Position',tp);
%-- yticks
YLIM_NTICKS = 5;
YLIM_SIG_FIGS_I = 2;
YLIM_SIG_FIGS = 2;
STD_FACTOR = 1;
%-- calculate ylim
if ~isempty(PLOT_STRUCT.ylim)
    y_lims = PLOT_STRUCT.ylim;
else
    y_lims = [min(table_in.(measure_char))-std(table_in.(measure_char))*STD_FACTOR,...
                max(table_in.(measure_char))+std(table_in.(measure_char))*STD_FACTOR];
    y_lims(2) = max([y_stats,y_lims(2)]); % stats correction
end
%-- calculate yticks
if ~isempty(PLOT_STRUCT.ytick)
    yticks = PLOT_STRUCT.ytick;
else
    if y_lims(2) > 0
        y_lims(2) = round(y_lims(2),YLIM_SIG_FIGS,'significant','TieBreaker','fromzero');
    else
        y_lims(2) = round(y_lims(2),YLIM_SIG_FIGS,'significant','TieBreaker','tozero');
    end
    if y_lims(1) > 0
        y_lims(1) = round(y_lims(1),YLIM_SIG_FIGS,'significant','TieBreaker','tozero');
    else
        y_lims(1) = round(y_lims(1),YLIM_SIG_FIGS,'significant','TieBreaker','fromzero');
    end
    u = y_lims(2);
    l = y_lims(1);
    bb = (u-l)/(YLIM_NTICKS-1);
    bbr = round(bb,YLIM_SIG_FIGS_I,'significant');
    dbb = abs(bbr-bb)*(YLIM_NTICKS-1);
    l = l-dbb/2;
    u = u+dbb/2;
    yticks = linspace(l,u,YLIM_NTICKS);
    y_lims = [min(yticks)-1e-8,max(yticks)+1e-8];
end
%-- ytick labels
if ~all(cellfun(@isempty,PLOT_STRUCT.ytick_labs))
    ytick_labs = PLOT_STRUCT.ytick_labs;
else
    tmp = yticks;
    tmp(tmp < 1e-8 & tmp > -1e-8) = 0;
    ytick_labs = cellstr(string(tmp));
end
%-- sets
set(ax,'YLim',y_lims, ...
    'YTick',yticks, ...
    'YTickLabel',ytick_labs);

%## AX TITLE
title(ax,PLOT_STRUCT.title, ...
    title_props{:});    

%## SET GROUP LABELS
axylim = get(ax,'YLim');
if PLOT_STRUCT.do_group_labels
    cnt_g = 1;
    for i = 1:length(groups)
        xx = gticks(i);
        yy = PLOT_STRUCT.group_label_offset(2)*abs(axylim(2)-axylim(1));
        tt = text(ax,xx,yy,0,char(PLOT_STRUCT.group_labels(i)),...
            group_label_props{:});
        tt.Units = 'normalized';
        tt.Position(2) = PLOT_STRUCT.group_label_offset(2);
        tt.Position(1) = tt.Position(1)+PLOT_STRUCT.group_label_offset(1);
        cnt_g = cnt_g + length(conds);
    end        
end
end