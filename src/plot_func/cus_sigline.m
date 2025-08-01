function [h_out,sig_gap] = cus_sigline(ax,p_value,conn_x,varargin)
%CUS_SIGLINE Summary of this function goes here
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
axylim = get(ax,'YLim');
%-
DEF_LINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_levels',[0.05,0.01,0.001],...
    'sig_offset_x',0,...
    'sig_offset_y',0); 
%## CHECK FUNCTIONS
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ax',@isgraphics);
addRequired(p,'p_value',@isnumeric);
addRequired(p,'conn_x',@isnumeric);
%- PARAMETER
addParameter(p,'LINE_STRUCT',DEF_LINE_STRUCT,@(x) validate_struct(x,DEF_LINE_STRUCT));
parse(p,ax,p_value,conn_x,varargin{:});
%- SET DEFAULTS
LINE_STRUCT = p.Results.LINE_STRUCT;
LINE_STRUCT = set_defaults_struct(LINE_STRUCT,DEF_LINE_STRUCT);
%% ===================================================================== %%
%- REDUNDANCY
LINE_STRUCT.sig_levels = sort(LINE_STRUCT.sig_levels,'descend');
%-- connector offset
if ~isempty(LINE_STRUCT.conn_y) && length(LINE_STRUCT.conn_y) == 2
    conn_y = LINE_STRUCT.conn_y;
else
    conn_y = repmat(gety(ax),[1,2]);
end
if ~isempty(LINE_STRUCT.conn_offset_y)
    conn_offset_y = LINE_STRUCT.conn_offset_y;
else
    conn_offset_y = 0.1*(abs(axylim(2)-axylim(1)));
end
conn_y = conn_y+conn_offset_y;
%-- signficance sign offset
if ~isempty(LINE_STRUCT.sig_offset_y)
    sig_offset_y = LINE_STRUCT.sig_offset_y;    
else
    sig_offset_y = abs(0.01*conn_y); %0.025*(abs(axylim(2)-axylim(1)));
end

%## PLOT
% sig_gap = [];
sig_gap = abs(sig_offset_y);
h_out = [];
hold on;
%- add label
chk1 = find((p_value <= LINE_STRUCT.sig_levels),1,'last');
%-- plot line
if ~isempty(chk1)    
    h_out = plot(ax,conn_x,conn_y+conn_offset_y,LINE_STRUCT.line_specs{:});
end
%-- add label
if ~isempty(chk1) && ~strcmp(LINE_STRUCT.sig_sign,'')    
    sig_sign = repmat(LINE_STRUCT.sig_sign,[1,chk1]);
    text(ax,mean(conn_x)+LINE_STRUCT.sig_offset_x,mean(conn_y)+sig_offset_y,sig_sign,...
        LINE_STRUCT.text_specs{:});
    % h_out = [h_out,tx];
end
hold on;
end

