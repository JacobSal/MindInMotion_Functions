function [ax,Pa,Li] = plot_psd(ax,psd_dat,freqs,varargin)
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
tt = tic();
%## DEFINE DEFAULTS
%- flags
do_split_colors = false;
do_print_chks = true;
Pa = [];
Li = [];

%## PLOT STRUCT
label_struct = struct('Units','normalized',...
    'FontName','Arial', ...
    'FontSize',12, ...
    'FontWeight','bold');
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
psd_line_struct = struct('LineWidth',2, ...
    'LineStyle','-', ...
    'DisplayName','line', ...
    'Color',[0.5,0.5,0.5,0.65] ...
    );
err_bnd_struct = struct( ...
    'LineStyle',':', ...
    'LineWidth',3, ...
    'FaceAlpha',0.6, ...
    'EdgeColor','none', ...
    'FaceColor',[0.5,0.5,0.5]);

%##
DEF_PLOT_STRUCT = struct( ...
    'do_set_ax_props',true, ...
    'ylim',[],...
    'y_label',{'10*log_{10}(PSD)'},...
    'y_label_props',label_struct, ...
    'y_label_offset',[0,0],...
    'ytick_labs',{{}}, ...
    'yticks',[], ...
    'xlim',[],...
    'x_label',{'Frequency (Hz)'},...
    'x_label_props',label_struct, ...
    'x_label_offset',[0,0],...
    'xtick_labs',{{}}, ...
    'xticks',[], ...
    'err_props',{{ ...
        'LineStyle',':', ...
        'LineWidth','3', ...
        'FaceAlpha',0.6, ...
        'EdgeColor','none', ...
        'FaceColor',[0.5,0.5,0.5]}}, ...
    'xtick_angle',45, ...
    'title',{{''}},...
    'title_props',title_props_struct, ...
    'title_offset',[0,0], ...
    'ax_props',ax_struct);
DEF_LINE_STRUCT = struct('do_line_avg',false, ...
    'line_props',psd_line_struct, ...
    'line_avg_fcn',@(x) mean(x,1), ...
    'do_err_shading',true, ...
    'err_bnd_props',err_bnd_struct, ...
    'err_bnd_vec',[]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ax',@(x) isgraphics(x))
addRequired(p,'psd_dat',@isnumeric);
addRequired(p,'freqs',@isnumeric);
%## PARAMETER
addParameter(p,'PLOT_STRUCT',DEF_PLOT_STRUCT,@(x) validate_struct(x,DEF_PLOT_STRUCT,do_print_chks));
addParameter(p,'LINE_STRUCT',DEF_LINE_STRUCT,@(x) validate_struct(x,DEF_LINE_STRUCT,do_print_chks));
%##
parse(p,ax,psd_dat,freqs,varargin{:});
%## SET DEFAULTS
PLOT_STRUCT = p.Results.PLOT_STRUCT;
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEF_PLOT_STRUCT,do_print_chks);
LINE_STRUCT = p.Results.LINE_STRUCT;
LINE_STRUCT = set_defaults_struct(LINE_STRUCT,DEF_LINE_STRUCT,do_print_chks);


%## LINE COLOR FOR MULTIPLE LINES
LINE_STRUCT.line_props.Color
if ~isempty(LINE_STRUCT.line_props.Color) && size(LINE_STRUCT.line_props.Color,1) > 1
    do_split_colors = true;
    color_hold = LINE_STRUCT.line_props.Color;
    LINE_STRUCT.line_props.Color = [];
end

%## CONVERT STRUCT ARGUMENTS TO CELL ARGS
line_props = struct2args(LINE_STRUCT.line_props);
y_label_props = struct2args(PLOT_STRUCT.y_label_props);
x_label_props = struct2args(PLOT_STRUCT.x_label_props);
title_props = struct2args(PLOT_STRUCT.title_props);
ax_props = struct2args(PLOT_STRUCT.ax_props);
%% ===================================================================== %%
hold on;
%## PLOT TRACE && ERROR SHADING
if LINE_STRUCT.do_err_shading
    tmp = LINE_STRUCT.err_bnd_vec;
    sz = size(tmp);
    if sz(2) ~= 2
        tmp = LINE_STRUCT.err_bnd_vec';
    end
    sz = size(psd_dat);
    if all(sz > 1) && LINE_STRUCT.do_line_avg
        [ax,Pa,Li] = cus_jackknife_sung(ax,freqs, ...
            LINE_STRUCT.line_avg_fcn(psd_dat), ...
            tmp, ...
            'LINE_STRUCT',LINE_STRUCT); %#ok<ASGLU>   
    else
         [ax,Pa,Li] = cus_jackknife_sung(ax,freqs, ...
            psd_dat, ...
            tmp, ...
            'LINE_STRUCT',LINE_STRUCT); %#ok<ASGLU>
    end
    hold on;
end

%## LINE
sz = size(psd_dat);
if all(sz > 1) && LINE_STRUCT.do_line_avg
    Li = plot(freqs,LINE_STRUCT.line_avg_fcn(psd_dat), ...
        line_props{:});
else
     Li = plot(freqs,psd_dat, ...
         line_props{:});
     if do_split_colors
         for i = 1:length(Li)
             set(Li(i),'Color',color_hold(i,:));
         end
     end
end

%## STATS
% [axsignif,Pa] = plot_psd_stats(ax, freqs, pcond, 'background', 'Frequency (Hz)');

%## set figure color, units, and size
% ylabel(ax,PLOT_STRUCT.y_label);
% xlabel(ax,PLOT_STRUCT.x_label,PLOT_STRUCT.x_label_props{:});
% legend(Li);
if PLOT_STRUCT.do_set_ax_props
    %-- set axes units, and size if they haven't already been
    if isempty(ax)
        ax = gca;
    end
    set(ax,ax_props{:});
    %-- xlabel
    xh = xlabel(ax,PLOT_STRUCT.x_label, ...
        x_label_props{:}); 
    %-- xticks
    if ~isempty(PLOT_STRUCT.xlim)
        set(ax,'XLim',PLOT_STRUCT.xlim)
    end
    if ~isempty(PLOT_STRUCT.xticks)
        set(ax,'XTick',sort(PLOT_STRUCT.xticks));
    end
    if ~isempty(PLOT_STRUCT.xtick_labs)
        set(ax,'XTickLabel',PLOT_STRUCT.xtick_labs);
    end
    xtickangle(PLOT_STRUCT.xtick_angle);
    %-- set pos
    tp = get(xh,'Position');
    tp(1,2) = tp(1,2)+PLOT_STRUCT.x_label_offset(2);
    tp(1,1) = tp(1,1)+PLOT_STRUCT.x_label_offset(1);
    set(xh,'Position',tp);    
    
    %-- set ylabel
    yh = ylabel(ax,PLOT_STRUCT.y_label, ...
        y_label_props{:});
    %-- xticks
    if ~isempty(PLOT_STRUCT.ylim)
        set(ax,'YLim',PLOT_STRUCT.ylim)
    end
    if ~isempty(PLOT_STRUCT.yticks)
        set(ax,'YTick',sort(PLOT_STRUCT.yticks));
    end
    if ~isempty(PLOT_STRUCT.ytick_labs)
        set(ax,'YTickLabel',PLOT_STRUCT.ytick_labs);
    end
    %-- set pos
    tp = get(yh,'Position');
    tp(1,2) = tp(1,2)+PLOT_STRUCT.y_label_offset(2);
    tp(1,1) = tp(1,1)+PLOT_STRUCT.y_label_offset(1);
    set(yh,'Position',tp);

    %-- title
    th = title(ax,PLOT_STRUCT.title, ...
        title_props{:});
    %-- set pos
    tp = get(th,'Position');
    tp(1,2) = tp(1,2)+PLOT_STRUCT.title_offset(2);
    tp(1,1) = tp(1,1)+PLOT_STRUCT.title_offset(1);
    set(th,'Position',tp);
    %-- time
end
fprintf('done. plot_psd.m %0.2g s\n',toc(tt));
end