function [ax] = plot_contourf_cell(ax,allersp,alltimes,allfreqs,...
    allersp_pcond,varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
cat_logo();
%-
PRINT_CHKS = true;
tmpfg = gcf();
psr = tmpfg.Position(4)/tmpfg.Position(3); % h/w (e.g., 9/6.5 for letter paper
%-
label_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
ax_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','normal', ...
    'OuterPosition',[0,0,1,1], ...
    'LineWidth',1, ...
    'Position',[0.06,0.7,0.13,0.16]);
DEFAULT_PLOT_STRUCT = struct( ...
    'title',{'ERSP'},...
    'title_props',label_struct, ...
    'stat_alpha',0.05, ...
    'xlabel','Gait Cycle Time (ms)',...
    'xlabel_props',label_struct, ...
    'ylabel','Frequency (Hz)',...
    'ylabel_props',label_struct, ...
    'yticks',[], ...
    'yscale',{'log'}, ...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'xticklabel_angle',45,...
    'ax_props',ax_struct, ...
    'clim',[],...
    'freq_lims',[],...
    'time_lims',[],...    
    'contourf_grain',30,...
    'alpha_multiple',0.6,...
    'do_display_sgtitle',true, ...
    'sgtitle_char',{''}, ...
    'sgtitle_shift',[0,0.65], ...
    'sgtitle_boxsz',[0.1,0.1], ...
    'sgtitle_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',12,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top',...
        'Units','normalized'), ...
    'do_display_colorbar',true, ...
    'colormap',linspecer(), ...
    'cbar_shift',[psr*0,psr*0,1,1],...
    'cbar_label_shift',[psr*0.2,(1/psr)*1.3,1,1],...
    'cbar_ticks',[],...
    'cbar_ticklabs',{{''}},...
    'cbar_label','\Delta Power (dB)',...
    'cbar_tickangle',0, ...
    'cbar_props',struct(...
        'Location','EastOutside', ...
        'FontSize',8, ...
        'FontWeight','bold'), ...
    'cbar_label_props',struct( ...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',8,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle',...
        'Units','normalized',...
        'Rotation',270), ...
    'do_display_bandmarks',true,...
    'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
    'bandmarks_shift',[0.075*psr,-0.05*(1/psr); ...
                       0.075*psr,0.25*(1/psr); ...
                       0.075*psr,0.53*(1/psr); ...
                       0.075*psr,0.89*(1/psr)],...
    'bandmarks_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',8,...
        'FontWeight','bold',...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Units','normalized'));

%% Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ax',@(x) isgraphics(x) || isempty(x));
addRequired(p,'allersp',@isnumeric);
addRequired(p,'alltimes',@isnumeric);
addRequired(p,'allfreqs',@isnumeric);
addRequired(p,'allersp_pcond',@isnumeric);
%## PARAMETER
addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT,PRINT_CHKS));
parse(p,ax,allersp,alltimes,allfreqs,allersp_pcond,varargin{:});
%## SET DEFAULTS
PLOT_STRUCT = p.Results.PLOT_STRUCT;
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT,PRINT_CHKS);

%## GET FREQUENCY & TIME LIMITS
if isempty(PLOT_STRUCT.freq_lims)
    PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
end
if isempty(PLOT_STRUCT.time_lims)
    PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
end

%## CONVERT STRUCT ARGUMENTS TO CELL ARGS
ylabel_props = struct2args(PLOT_STRUCT.ylabel_props);
xlabel_props = struct2args(PLOT_STRUCT.xlabel_props);
title_props = struct2args(PLOT_STRUCT.title_props);
bandmarks_props = struct2args(PLOT_STRUCT.bandmarks_props);
cbar_label_props = struct2args(PLOT_STRUCT.cbar_label_props);
sgtitle_props = struct2args(PLOT_STRUCT.sgtitle_props);
ax_props = struct2args(PLOT_STRUCT.ax_props);
cbar_props = struct2args(PLOT_STRUCT.cbar_props);
%% COLORLIMITS ALG
%## set color limits
if isempty(PLOT_STRUCT.clim)
    %##
    data = allersp(:);
    % bound = max([abs(prctile(data,5,'all')),abs(prctile(data,95,'all'))]);
    bound = [prctile(data,5,'all'),prctile(data,95,'all')];
    %-- assign limit and c-lim labels
    fudge = bound*0.1;
    % PLOT_STRUCT.clim = sort([-round(bound,1,'significant'),round(bound,1,'significant')]);
    PLOT_STRUCT.clim = sort([round(bound(1),1,'significant'),round(bound(2),1,'significant')]);
    if length(unique(PLOT_STRUCT.clim)) == 1
        PLOT_STRUCT.clim = [PLOT_STRUCT.clim(1)-fudge,PLOT_STRUCT.clim(2)+fudge];
    end
    PLOT_STRUCT.cbar_ticks = round(linspace(PLOT_STRUCT.clim(1),PLOT_STRUCT.clim(2),5),2,'significant');
    PLOT_STRUCT.cbar_ticklabs = cellstr(string(PLOT_STRUCT.cbar_ticks));
end
%%
inds1 = allfreqs>=PLOT_STRUCT.freq_lims(1) & allfreqs<=PLOT_STRUCT.freq_lims(2);
inds2 = alltimes>=PLOT_STRUCT.time_lims(1) & alltimes<=PLOT_STRUCT.time_lims(2);
tmp_freqs = allfreqs(inds1);
tmp_times = alltimes(inds2);
if size(allersp,1) ~= length(tmp_freqs) || size(allersp,2) ~= length(tmp_times)
    allersp = allersp(inds1,inds2,:);
    allersp_pcond = allersp_pcond(inds1,inds2,:);
    allfreqs = tmp_freqs;
    alltimes = tmp_times;
end
PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];

%%
%## AXES PLOT
if isempty(ax)
    ax = axes();
end
hold on;
%-- mask all non-significant content
allersp_bin = allersp_pcond > PLOT_STRUCT.stat_alpha;
% allersp_mask = allersp.*allersp_bin;
alpha_mask = ones(size(allersp_pcond))*PLOT_STRUCT.alpha_multiple;
alpha_mask(~allersp_bin) = 0;
alpha_mask(allersp_bin) = PLOT_STRUCT.alpha_multiple;
%--
contourf(alltimes,allfreqs,allersp,PLOT_STRUCT.contourf_grain,...
           'linecolor','none');
hold on;
% allersp_mask(allersp_mask == 0) = -1e10;
% imagesc(alltimes,allfreqs,allersp_mask,'AlphaData',alpha_mask);
cmap_ct = gray();
c_stat_val = [1,1,1];
for i = 1:size(alpha_mask,3)
    imin = allersp_bin(:,:,i)*c_stat_val(i);
    scaledpict = round(1 + mat2gray(imin)*(size(cmap_ct,1)-1));
    imin = uint8(ind2rgb(scaledpict,cmap_ct));
    image(alltimes,allfreqs,imin, ...
        'AlphaData',alpha_mask(:,:,i)*0.7);
end
% if size(alpha_mask,3) > 1
%     % cmap_ct = linspecer();
%     cmap_ct = gray();
%     c_stat_val = [1,0.3,0.5];
%     for i = 1:size(alpha_mask,3)
%         imin = allersp_bin(:,:,i)*c_stat_val(i);
%         scaledpict = round(1 + mat2gray(imin)*(size(cmap_ct,1)-1));
%         imin = uint8(ind2rgb(scaledpict,cmap_ct));
%         image(alltimes,allfreqs,imin, ...
%             'AlphaData',alpha_mask(:,:,i)*0.7);
%     end
% else
%     image(alltimes,allfreqs,uint8(repmat(allersp_bin,1,1,3)), ...
%         'AlphaData',alpha_mask);
% end
colormap(PLOT_STRUCT.colormap);
%-
set(ax,ax_props{:}); %'LineWidth',1, ...
    % 'FontName',PLOT_STRUCT.font_name,...
    % 'FontSize',PLOT_STRUCT.font_size,...
    % 'FontWeight','bold',...
    % 'OuterPosition',[0 0 1 1]);
set(ax,'CLim',PLOT_STRUCT.clim,...
    'XLim',PLOT_STRUCT.time_lims,...
    'YDir','norm',...
    'YLim',PLOT_STRUCT.freq_lims,...
    'YScale',PLOT_STRUCT.yscale, ...
    'XScale','log')


%## SET YTICKS & YLABELS
yticks = [];
ytick_labels = [];
MAX_ERR = 0.95;
if isempty(PLOT_STRUCT.yticks)
    if PLOT_STRUCT.freq_lims(2) <= 80
        chk = abs(PLOT_STRUCT.freq_lims(2)-80) < MAX_ERR;
        if chk
            mv = PLOT_STRUCT.freq_lims(2);
        else
            mv = 80;
        end
        yticks = [4.01,8,13,30,50,mv];
        ytick_labels = {'4','8','13','30','50','80'};
    elseif PLOT_STRUCT.freq_lims(2) <= 100
        yticks = [4.01,8,13,30,50,99.4843];
        ytick_labels = {'4','8','13','30','50','100'};    
    elseif PLOT_STRUCT.freq_lims(2) <= 200
        yticks = [4.01,8,13,30,50,99.4843,195.5];
        ytick_labels = {'4','8','13','30','50','100','200'};
    end
    if ~isempty(yticks) && ~isempty(ytick_labels)
        set(ax,'YTick',yticks, ...
            'YTickLabel',ytick_labels);
        for k = 1:length(yticks)
            yline(ax,yticks(k),'k:');
        end
    end
else
    set(ax,'YTick',PLOT_STRUCT.yticks)
end
%-- label
ylabel(PLOT_STRUCT.ylabel, ...
    ylabel_props{:});

%## SET XTICKS & XLABELS
%- set x-axis ticks
if ~isempty(PLOT_STRUCT.xticklabel_times)
    xrng = get(ax,'XLim');
    if PLOT_STRUCT.xticklabel_times(1) < xrng(1)
        PLOT_STRUCT.xticklabel_times(1) = xrng(1);
    end
    if PLOT_STRUCT.xticklabel_times(end) > xrng(end)
        PLOT_STRUCT.xticklabel_times(end) = xrng(end);
    end
    set(ax,'XTick',PLOT_STRUCT.xticklabel_times);
    for k = 1:length(PLOT_STRUCT.xticklabel_times)
        xline(ax,PLOT_STRUCT.xticklabel_times(k),'k--')
    end
end
if ~isempty(PLOT_STRUCT.xticklabel_chars)
    set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
    xtickangle(PLOT_STRUCT.xticklabel_angle);
end
%-- label
xlabel(PLOT_STRUCT.xlabel, ...
    xlabel_props{:});

%## SET AXES TITLE
title(PLOT_STRUCT.title, ...
    title_props{:});

%- set color bar
set(ax,'Position',PLOT_STRUCT.ax_props.Position);
ANN_DIM = repmat(min([PLOT_STRUCT.ax_props.Position(3),PLOT_STRUCT.ax_props.Position(4)]),1,2);
% ANN_DIM = [0.1,0.1];
% npos = get(ax,'Position');
if PLOT_STRUCT.do_display_colorbar    
    c = colorbar(ax);
    set(c,cbar_props{:});
    xx = c.Position(1)+PLOT_STRUCT.cbar_shift(1)*PLOT_STRUCT.ax_props.Position(3);
    yy = c.Position(2)+PLOT_STRUCT.cbar_shift(2)*PLOT_STRUCT.ax_props.Position(4);
    c.Ruler.TickLabelRotation = PLOT_STRUCT.cbar_tickangle;
    % c.Limits = PLOT_STRUCT.clim;
    if isempty(PLOT_STRUCT.cbar_ticks) || isempty(PLOT_STRUCT.cbar_ticklabs)
        PLOT_STRUCT.cbar_ticklabs = linspace(PLOT_STRUCT.clim(1),PLOT_STRUCT.clim(2),5);
        PLOT_STRUCT.cbar_ticks = linspace(PLOT_STRUCT.clim(1),PLOT_STRUCT.clim(2),5);
    end
    set(c,'Position',[xx,yy,c.Position(3)*PLOT_STRUCT.cbar_shift(3),c.Position(4)*PLOT_STRUCT.cbar_shift(4)], ...
        'Limits',PLOT_STRUCT.clim, ...
        'Ticks',PLOT_STRUCT.cbar_ticks, ...
        'TickLabels',PLOT_STRUCT.cbar_ticklabs);
    %- color bar label
    % ANN_DIM = [c.Position(3),c.Position(4)];
    % ANN_DIM = repmat(min([PLOT_STRUCT.ax_props.Position(3),PLOT_STRUCT.ax_props.Position(4)]),1,2);
    xx = c.Position(1)+ANN_DIM(1)*PLOT_STRUCT.cbar_label_shift(1);
    yy = c.Position(2)+ANN_DIM(2)*PLOT_STRUCT.cbar_label_shift(2);
    % xx = c.Position(1)*PLOT_STRUCT.cbar_label_shift(1)+(ANN_DIM(1)/2);
    % yy = c.Position(2)*PLOT_STRUCT.cbar_label_shift(2)+(ANN_DIM(2)/2);
    
    % disp([xx,yy])
    a = annotation(gcf,'textbox',...        
        [xx,yy,ANN_DIM.*PLOT_STRUCT.cbar_label_shift(3:4)],...
        'String',PLOT_STRUCT.cbar_label, ...
        cbar_label_props{:});
end

if PLOT_STRUCT.do_display_bandmarks
    % ANN_DIM = [0.1,0.1];
    npos = get(ax,'Position');
    for ii = 1:length(PLOT_STRUCT.bandmarks)
        xx = npos(1)+npos(3)*PLOT_STRUCT.bandmarks_shift(ii,1)-ANN_DIM(1)/2;
        yy = npos(2)+npos(4)*PLOT_STRUCT.bandmarks_shift(ii,2)-ANN_DIM(2)/2;
        a = annotation(gcf,'textbox',[xx,yy,ANN_DIM],...
            'String',PLOT_STRUCT.bandmarks{ii}, ...
            bandmarks_props{:});
    end
end

if PLOT_STRUCT.do_display_sgtitle
    xx = 0.5+(-PLOT_STRUCT.sgtitle_boxsz(1)/2)+PLOT_STRUCT.sgtitle_shift(1); % always aim for center of figure
    yy = PLOT_STRUCT.ax_props.Position(2)+PLOT_STRUCT.ax_props.Position(4)*PLOT_STRUCT.sgtitle_shift(2);
    a = annotation(gcf,'textbox',[xx,yy,PLOT_STRUCT.sgtitle_boxsz],...
        'String',PLOT_STRUCT.sgtitle_char, ...
        sgtitle_props{:});
end
set(ax,'Position',PLOT_STRUCT.ax_props.Position);
end
%--
function [argout,ind] = get_cellarg(cellargs,argchar)
    ind = find(strcmp(argchar,cellargs));
    argout = cellargs(ind+1);
end