function [fig_handles] = imu_valid_plots(EEG,varargin)
% function [fig_handles] = imu_valid_plots(data_in,data_labs,varargin)
%IMU_VALID_PLOTS Summary of this function goes here
%   IN: 
%       data_in, DOUBLE
%           matrix of size (NxMxK) where N is the number of data channels,
%           M is the number data points, and K is the number of groups
%           (e.g., body positions, world positions, sensor positions)
%       data_labs, CELL
%           cell of chars of size(NxK) where N is the labels the channel
%           for every N in group K. 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%## GENERALIZED SIGNAL PLOTTING FOR IMU
% %- logo
% cat_logo();
% %## DEFAULTS
% AXES_DEFAULT_PROPS = {'box','off', ...
%     'xtick',[], ...
%     'ytick',[], ...
%     'ztick',[], ...
%     'xcolor',[1,1,1], ...
%     'ycolor',[1,1,1]};
% save_fname_ext = [];
% save_dir = [];
% %## PARSE
% p = inputParser;
% %## REQUIRED
% addRequired(p,'data_in',@isnumeric);
% addRequired(p,'data_labs',@iscell);
% addOptional(p,'save_dir',@isstruct);
% addRequired(p,'save_fname_ext',@isstruct);
% %## OPTIONAL
% %## PARAMETER
% parse(p, data_in, varargin{:});
% %## SET DEFAULTS
% save_dir = p.Results.save_dir;
% save_fname_ext = p.Results.save_fname_ext;
% fig_handles = [];

%## QUICK/DIRTY PLOT IMPLEMENTATION TO BE PAIRED WITH THE
%"IMU_GET_POS_COORDS.M" FUNCTION
%- logo
cat_logo();
%## DEFAULTS
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[], ...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
save_fname_ext = [];
save_dir = [];
EXPORT_RES = 300;
DO_LIMIT_AXES = false;
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
addOptional(p,'save_dir',save_dir,@(x) ischar(x) || isempty(x));
addOptional(p,'save_fname_ext',save_fname_ext,@(x) ischar(x) || isempty(x));
%## PARAMETER
addParameter(p,'EXPORT_RES',EXPORT_RES,@isnumeric);
addParameter(p,'DO_LIMIT_AXES',DO_LIMIT_AXES,@isnumeric);
parse(p, EEG, varargin{:});
%## SET DEFAULTS
save_dir = p.Results.save_dir;
save_fname_ext = p.Results.save_fname_ext;
EXPORT_RES = p.Results.EXPORT_RES;
DO_LIMIT_AXES = p.Results.DO_LIMIT_AXES;
fig_handles = [];
%% ===================================================================== %%
%## EXTRACT DATA
%-
body_posx = find(contains({EEG.chanlocs.labels},'body_xpos','IgnoreCase',true));
body_posy = find(contains({EEG.chanlocs.labels},'body_ypos','IgnoreCase',true));
body_posz = find(contains({EEG.chanlocs.labels},'body_zpos','IgnoreCase',true));
world_posx = find(contains({EEG.chanlocs.labels},'world_xpos','IgnoreCase',true));
world_posy = find(contains({EEG.chanlocs.labels},'world_ypos','IgnoreCase',true));
world_posz = find(contains({EEG.chanlocs.labels},'world_zpos','IgnoreCase',true));
sensor_posx = find(contains({EEG.chanlocs.labels},'sensor_xpos','IgnoreCase',true));
sensor_posy = find(contains({EEG.chanlocs.labels},'sensor_ypos','IgnoreCase',true));
sensor_posz = find(contains({EEG.chanlocs.labels},'sensor_zpos','IgnoreCase',true));

%- check continuous or epoch
if length(size(EEG.data)) == 3
    span = 10; % grab 5 strides to display
    s_i = randi(size(EEG.data,3)-(span+1),1);
    s_i = s_i:s_i+span;
    imu_times = EEG.times;
    %- chan data
    body_posx = squeeze(EEG.data(body_posx,:,s_i));
    body_posy = squeeze(EEG.data(body_posy,:,s_i));
    body_posz = squeeze(EEG.data(body_posz,:,s_i));
    world_posx = squeeze(EEG.data(world_posx,:,s_i));
    world_posy = squeeze(EEG.data(world_posy,:,s_i));
    world_posz = squeeze(EEG.data(world_posz,:,s_i));
    sensor_posx = squeeze(EEG.data(sensor_posx,:,s_i));
    sensor_posy = squeeze(EEG.data(sensor_posy,:,s_i));
    sensor_posz = squeeze(EEG.data(sensor_posz,:,s_i));
    %-
    cmaps = linspecer(length(s_i));
    is_cont = false;
elseif length(size(EEG.data)) == 2
    if size(EEG.data) < 100000
        s_i = 1:size(EEG.data,2);
        imu_times = EEG.times(s_i);
    else
        span = 10000; % display 10000 ms of data
        s_i = randi(size(EEG.data,2)-(span+1),1);
        s_i = s_i:s_i+span;
        imu_times = EEG.times(s_i);
    end
    %- chan data
    body_posx = squeeze(EEG.data(body_posx,s_i));
    body_posy = squeeze(EEG.data(body_posy,s_i));
    body_posz = squeeze(EEG.data(body_posz,s_i));
    world_posx = squeeze(EEG.data(world_posx,s_i));
    world_posy = squeeze(EEG.data(world_posy,s_i));
    world_posz = squeeze(EEG.data(world_posz,s_i));
    sensor_posx = squeeze(EEG.data(sensor_posx,s_i));
    sensor_posy = squeeze(EEG.data(sensor_posy,s_i));
    sensor_posz = squeeze(EEG.data(sensor_posz,s_i));
    cmaps = linspecer(3);
    is_cont = true;
else
    error('Size of data is not understood for this function...')
end
%- reflect y position to be anatomicallly appropriate?
% body_posy = -body_posy;
%- get ground reaction force
grfl = find(contains({EEG.chanlocs.labels},'GRF_L','IgnoreCase',true));
grfr = find(contains({EEG.chanlocs.labels},'GRF_R','IgnoreCase',true));
%-
if ~isempty(grfl) && ~isempty(grfr)
    if is_cont
        grfl = squeeze(EEG.data(grfl,s_i));
        grfr = squeeze(EEG.data(grfr,s_i));
    else
        grfl = squeeze(EEG.data(grfl,:,s_i));
        grfr = squeeze(EEG.data(grfr,:,s_i));
    end
end
%% PLOT 1
pp1 = [];
pp2 = [];
pp3 = [];
line_alpha = 0.75;
%-
LINE_PROPS3D_BODY = {'LineWidth',1, ...
    'LineStyle','-', ...
    'DisplayName','body frame'};
LINE_PROPS3D_WORLD = {'LineWidth',1, ...
    'LineStyle','-', ...
    'DisplayName','world frame'};
LINE_PROPS3D_SENSOR = {'LineWidth',1, ...
    'LineStyle','-', ...
    'DisplayName','sensor frame'};
%-
fig = figure('color','white');
set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
set(gca,AXES_DEFAULT_PROPS{:})
ax = axes();
hold on;
if is_cont
    % pp1 = plot3(ax,sensor_posx,sensor_posy,sensor_posz,LINE_PROPS3D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot3(ax,body_posx,body_posy,body_posz,LINE_PROPS3D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot3(ax,world_posx,world_posy,world_posz,LINE_PROPS3D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posx,2)
        hold on;
        % pp1 = plot3(ax,sensor_posx(:,s_i),sensor_posy(:,s_i),sensor_posz(:,s_i),LINE_PROPS3D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot3(ax,body_posx(:,s_i),body_posy(:,s_i),body_posz(:,s_i),LINE_PROPS3D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot3(ax,world_posx(:,s_i),world_posy(:,s_i),world_posz(:,s_i),LINE_PROPS3D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
if DO_LIMIT_AXES
    xlim([min(body_posx)-0.5*std(body_posx),max(body_posx)+0.5*std(body_posx)])
    ylim([min(body_posy)-0.5*std(body_posy),max(body_posy)+0.5*std(body_posy)])
    zlim([min(body_posz)-0.5*std(body_posz),max(body_posz)+0.5*std(body_posz)])
end
if isempty(pp1)
    legend([pp2,pp3])
else
    legend([pp1,pp2,pp3])
end
view(ax,[45,25])
xlabel('Anteroposterior (m)');
ylabel('Mediolateral (m)');
zlabel('Up-down (m)');
set(ax,'Position',[0.25,0.5,0.5,0.5]);  %[left bottom width height]                
set(ax,'projection', 'perspective', 'box', 'on')
align_axislabel([], ax)

hold off;
exportgraphics(fig,[save_dir filesep sprintf('%s_3dstate_imu_fig.png',save_fname_ext)],'Resolution',EXPORT_RES)
%-
fig_handles = [fig_handles,fig];
%% PLOT 2
y_shift = 0.225;
ax_y = 0.75;
ax_x = 0.25;
ax_h = 0.175;
ax_w = 0.5;
leg_token_size = 10;
LINE_PROPS2D_GRFL = {'LineWidth',1, ...
    'LineStyle','-.', ...
    'DisplayName','GRF L'};
LINE_PROPS2D_GRFR = {'LineWidth',1, ...
    'LineStyle','-.', ...
    'DisplayName','GRF R'};
LINE_PROPS2D_BODY = {'LineWidth',1, ...
    'LineStyle','-.', ...
    'DisplayName','body frame'};
LINE_PROPS2D_WORLD = {'LineWidth',1, ...
    'LineStyle',':', ...
    'DisplayName','world frame'};
LINE_PROPS2D_SENSOR = {'LineWidth',1, ...
    'LineStyle',':', ...
    'DisplayName','sensor frame'};
%-
fig = figure('color','white');
set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
set(gca,AXES_DEFAULT_PROPS{:})
%-
if ~isempty(grfl) && ~isempty(grfr)
    ax = axes();
    tmpc = linspecer(2);
    hold on;
    if is_cont
        pp = plot(ax,imu_times,grfl,LINE_PROPS2D_GRFL{:});
        pp.Color = [tmpc(1,:),0.75];
        pp = plot(ax,imu_times,grfr,LINE_PROPS2D_GRFR{:});
        pp.Color = [tmpc(2,:),0.75];
    else
        for s_i = 1:size(body_posx,2)
            pp = plot(ax,imu_times,grfl(:,s_i),LINE_PROPS2D_GRFL{:});
            pp.Color = [tmpc(1,:),0.75];
            pp = plot(ax,imu_times,grfr(:,s_i),LINE_PROPS2D_GRFR{:});
            pp.Color = [tmpc(2,:),0.75];
        end
    end
    set(ax,'Position',[ax_x,ax_y,ax_w,ax_h]);  %[left bottom width height]
    ylabel(ax,'Force (N)');
    xlabel(ax,'Time (ms)');
end
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,imu_times,sensor_posy,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,imu_times,body_posy,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,imu_times,world_posy,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posy,2)
        % pp1 = plot(ax,imu_times,sensor_posy(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,imu_times,body_posy(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,imu_times,world_posy(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
% legend([pp1,pp2,pp3]);
set(ax,'Position',[ax_x,ax_y-y_shift,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Mediolateral Position (m)');
xlabel(ax,'Time (ms)');
% if DO_LIMIT_AXES
%     ylim([min(world_posy)-0.5*std(world_posy),max(world_posy)+0.5*std(world_posy)])
% end
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,imu_times,sensor_posx,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,imu_times,body_posx,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,imu_times,world_posx,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posx,2)
        % pp1 = plot(ax,imu_times,sensor_posx(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,imu_times,body_posx(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,imu_times,world_posx(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
set(ax,'Position',[ax_x,ax_y-y_shift*2,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Anteroposterior Position (m)');
xlabel(ax,'time (ms)');
% if DO_LIMIT_AXES
%     ylim([min(world_posx)-0.5*std(world_posx),max(world_posx)+0.5*std(world_posx)])
% end
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,imu_times,sensor_posz,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,imu_times,body_posz,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,imu_times,world_posz,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posz,2)
        % pp1 = plot(ax,imu_times,sensor_posz(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,imu_times,body_posz(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,imu_times,world_posz(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
set(ax,'Position',[ax_x,ax_y-y_shift*3,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Up-down Position (m)');
xlabel(ax,'time (ms)');
% if DO_LIMIT_AXES
%     ylim([min(world_posz)-0.5*std(world_posz),max(world_posz)+0.5*std(world_posz)])
% end
%-
if isempty(pp1)
    legend(ax,[pp2,pp3])
else
    legend(ax,[pp1,pp2,pp3])
end
[leg,~,~,~]  = legend('boxoff');
tmp = get(leg,'String');
set(leg,'String',tmp,'FontName','Arial','FontSize',10)
set(leg,'Orientation','horizontal')
set(leg,'Units','normalized')
set(leg,'Position',[ax_x-0.1,...
    ax_y-(y_shift*3+0.06),leg.Position(3),leg.Position(4)]);
leg.ItemTokenSize(1) = leg_token_size;
hold off;
exportgraphics(fig,[save_dir filesep sprintf('%s_2dtime_imu_fig.png',save_fname_ext)],'Resolution',EXPORT_RES)
%-
fig_handles = [fig_handles,fig];
%% PLOT 3
x_shift = 0.275;
ax_y = 0.25;
ax_x = 0.1;
ax_h = 0.5;
ax_w = 0.2;
leg_token_size = 10;
% LINE_PROPS2D_BODY = {'LineWidth',1, ...
%     'LineStyle','-.', ...
%     'DisplayName','body frame'};
% LINE_PROPS2D_WORLD = {'LineWidth',1, ...
%     'LineStyle',':', ...
%     'DisplayName','world frame'};
%-
fig = figure('color','white');
set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
set(gca,AXES_DEFAULT_PROPS{:})
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,sensor_posy,sensor_posx,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,body_posy,body_posx,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,world_posy,world_posx,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posz,2)
        % pp1 = plot(ax,sensor_posy,sensor_posx(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,body_posy(:,s_i),body_posx(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,world_posy(:,s_i),world_posx(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
set(ax,'Position',[ax_x,ax_y,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Anteroposterior Pos. (m)');
xlabel(ax,'Mediolateral Pos. (m)');
title(ax,'Transverse')
if DO_LIMIT_AXES
    xlim([min(body_posy)-0.5*std(body_posy),max(body_posy)+0.5*std(body_posy)])
    ylim([min(body_posx)-0.5*std(body_posx),max(body_posx)+0.5*std(body_posx)])
end
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,sensor_posx,sensor_posz,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,body_posx,body_posz,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,world_posx,world_posz,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posz,2)
        % pp1 = plot(ax,sensor_posx,sensor_posz(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,body_posx(:,s_i),body_posz(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,world_posx(:,s_i),world_posz(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
set(ax,'Position',[ax_x+x_shift,ax_y,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Up-down Pos. (m)');
xlabel(ax,'Anteroposterior Pos. (m)');
title(ax,'Sagittal');
if DO_LIMIT_AXES
    xlim(ax,[min(body_posx)-0.5*std(body_posx),max(body_posx)+0.5*std(body_posx)])
    ylim(ax,[min(body_posz)-0.5*std(body_posz),max(body_posz)+0.5*std(body_posz)])
end
%-
ax = axes();
hold on;
if is_cont
    % pp1 = plot(ax,sensor_posy,sensor_posz,LINE_PROPS2D_SENSOR{:});
    % pp1.Color = [cmaps(1,:),line_alpha];
    pp2 = plot(ax,body_posy,body_posz,LINE_PROPS2D_BODY{:});
    pp2.Color = [cmaps(2,:),line_alpha];
    pp3 = plot(ax,world_posy,world_posz,LINE_PROPS2D_WORLD{:});
    pp3.Color = [cmaps(3,:),line_alpha];
else
    for s_i = 1:size(body_posz,2)
        % pp1 = plot(ax,sensor_posy,sensor_posz(:,s_i),LINE_PROPS2D_SENSOR{:});
        % pp1.Color = [cmaps(s_i,:),line_alpha];
        pp2 = plot(ax,body_posy(:,s_i),body_posz(:,s_i),LINE_PROPS2D_BODY{:});
        pp2.Color = [cmaps(s_i,:),line_alpha];
        pp3 = plot(ax,world_posy(:,s_i),world_posz(:,s_i),LINE_PROPS2D_WORLD{:});
        pp3.Color = [cmaps(s_i,:),line_alpha];
    end
end
set(ax,'Position',[ax_x+x_shift*2,ax_y,ax_w,ax_h]);  %[left bottom width height]
ylabel(ax,'Up-down Position (m)');
xlabel(ax,'Mediolateral Pos. (m)');
title(ax,'Coronal');
if DO_LIMIT_AXES
    xlim(ax,[min(body_posy)-0.5*std(body_posy),max(body_posy)+0.5*std(body_posy)])
    ylim(ax,[min(body_posz)-0.5*std(body_posz),max(body_posz)+0.5*std(body_posz)])
end
%-
if isempty(pp1)
    legend(ax,[pp2,pp3])
else
    legend(ax,[pp1,pp2,pp3])
end
[leg,~,~,~]  = legend('boxoff');
tmp = get(leg,'String');
set(leg,'String',tmp,'FontName','Arial','FontSize',10)
set(leg,'Orientation','horizontal')
set(leg,'Units','normalized')
set(leg,'Position',[ax_x-0.1,...
    ax_y-(0.075),leg.Position(3),leg.Position(4)]);
leg.ItemTokenSize(1) = leg_token_size;

%-
hold off;
exportgraphics(fig,[save_dir filesep sprintf('%s_2dstate_imu_fig.png',save_fname_ext)],'Resolution',EXPORT_RES)
%-
fig_handles = [fig_handles,fig];
end

