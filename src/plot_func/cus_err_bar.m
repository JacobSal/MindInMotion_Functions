function [pp] = cus_err_bar(ax,x_val,y_bnds,varargin)
%CUS_SIGBRACKET Summary of this function goes here
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
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
axylim = get(ax,'YLim');
%-
DEF_ERRBAR_STRUCT = struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
    'err_bar_width',0.1);
%## PARSER
p = inputParser;
%- REQUIRED
addRequired(p,'ax',@isgraphics)
addRequired(p,'x_val',@isnumeric)
addRequired(p,'y_bnds',@isnumeric)
%- PARAMETER
addParameter(p,'ERRBAR_STRUCT',DEF_ERRBAR_STRUCT,@(x) validate_struct(x,DEF_ERRBAR_STRUCT));
parse(p,ax,x_val,y_bnds,varargin{:});
%- SET DEFAULTS
ERRBAR_STRUCT = p.Results.ERRBAR_STRUCT;
ERRBAR_STRUCT = set_defaults_struct(ERRBAR_STRUCT,DEF_ERRBAR_STRUCT);
%% ===================================================================== %%
%## CALC LINE VALUES
upr_hline = [y_bnds(2),y_bnds(2)];
lwr_hline = [y_bnds(1),y_bnds(1)];
x_hline = [x_val-ERRBAR_STRUCT.err_bar_width,x_val+ERRBAR_STRUCT.err_bar_width];
%--
vline = [y_bnds(1),y_bnds(2)];
x_vline = [x_val,x_val];
%## PLOT
hold on;
p1 = plot(ax,x_hline,upr_hline,ERRBAR_STRUCT.line_specs{:});
p2 = plot(ax,x_vline,vline,ERRBAR_STRUCT.line_specs{:});
p3 = plot(ax,x_hline,lwr_hline,ERRBAR_STRUCT.line_specs{:});
%-- patch implementation?
pp = [p1,p2,p3];
end