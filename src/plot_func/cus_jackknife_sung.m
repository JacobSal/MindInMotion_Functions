function [ax,Pa,Li] = cus_jackknife_sung(ax,xx,yy,err_bounds,varargin)
% JACKKNIFE plots jackknife errorbars around a given curve
%
%     [Pa,Li,t] = JACKKNIFE(x,y,L,U,'r','g')
%     Plots a gray Jackknife around the line displayed in black very useful for
%     funky nature style error bars which are shaded.
%
%         Pa is a patch object for more help on patch objects see below
%         Li is a line object, more help on line object is available in MATLAB
%
%     USAGE :
%              1)   [Pa,Li] = JackKnife(x,y,E)
%                    Calculates the Lower and upper errorbars as
%                    L = Y-E and U = Y+E. It then takes a default gray color
%                    as patch color, and the line color as black and plots
%                    it around the line using a patch object.
%
%              2)   [Pa,Li] = JackKnife(x,y,E,LineColor,PatchColor)
%                    Calculates the Lower and upper errorbars as
%                    L = Y-E and U = Y+E. It then takes PatchColor
%                    as patch color, and the Line Color from the LineColor
%                    variable. It then plots it around the line
%                    using a patch object.
%
%              3)   [Pa,Li] = JackKnife(x,y,L,U)
%                    User Supplied bounds are taken as L and U, It then takes
%                    a default gray color as patch color, and the line color
%                    as black and plots it around the line using a patch object.
%
%              4)   [Pa,Li] = JackKnife(x,y,L,U,LineColor,PatchColor)
%                    User Supplied bounds are taken as L and U, It then takes
%                    PatchColor as patch color, and the Line Color from the LineColor
%                    variable. It then plots it around the line using a
%                    patch object.
%      CAVEATS
%                 1) Can be Slow sometimes for length(Array) > 10000,
%                 2) Needs better vectorization
%      EXAMPLE
%                         t = [-5:0.05:5];
%                         Y = sin(t);
%                         E = 0.4*rand(1,length(t));
%                         [Pa,Li] = JackKnife(t,Y,E);
%                         xlabel('time');
%                         ylabel('Amplitude');
%                         title('Using Errors alone');
%
%                         figure;
%                         L = Y - E;
%                         U = Y + E;
%                         [Pa,Li] = JackKnife(t,Y,L,U);
%                         xlabel('time');
%                         ylabel('Amplitude');
%                         title('Using Lower and Upper Confidence Intervals');
%                         hold on;
%
%                         Y1 = 2*Y;
%                         L = Y1 - 0.2;
%                         U = Y1 + 0.2;
%                         [Pa,Li] = JackKnife(t,Y1,L,U,[255 51 51]./255,[255 153 102]./255);
%                         hold on;
%                         [Pa,Li] = JackKnife(t,Y1*2,E,[51 51 153]./255,[102 153 204]./255);
%                         [Pa,Li] = JackKnife(t,Y1*2,E,'r','g');
% See also ERRORBAR, PATCH, LINE
%
%
% Version 0.001 Chandramouli Chandrasekaran (Chandt) - 13 April 2006.
% Modified 04/12 by James Finley to automatically compute means and
% standard errors for matrices, cell arrays, and data structures.
%
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
DEF_LINE_STRUCT = struct('line_width',2, ...
    'line_style','-', ...
    'line_alpha',1, ...
    'line_color',[1,1,1], ...
    'line_label',{'label'}, ...
    'line_avg_fcn',@mean, ...
    'line_avg_fcn_params',{{2}}, ...
    'do_err_shading',true, ...
    'err_alpha',0.6, ...
    'err_color',[0.5,0.5,0.5], ...
    'err_edge_color',[], ...
    'err_upr_bnd_fcn',@(x,p) 1*std(x,p{:}), ...
    'err_lwr_bnd_fcn',@(x,p) -1*std(x,p{:}), ...
    'err_upr_bnd_params',{{[],2}}, ...
    'err_lwr_bnd_params',{{[],2}}, ...
    'err_line_style',':', ...
    'err_line_width',3);
%## CHECK FUNCTIONS
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ax',@isgraphics);
addRequired(p,'xx',@isnumeric);
addRequired(p,'yy',@isnumeric);
addRequired(p,'err_bounds',@(x) (isnumeric(x) && max(size(x))==length(xx)));
%- PARAMETER
addParameter(p,'LINE_STRUCT',DEF_LINE_STRUCT,@(x) validate_struct(x,DEF_LINE_STRUCT));
% addParameter(p,'PLOT_STRUCT',PLOT_STRUCT,@(x) validate_struct(x,DEF_PLOT_STRUCT));
parse(p,ax,xx,yy,err_bounds,varargin{:});
%- SET DEFAULTS
LINE_STRUCT = p.Results.LINE_STRUCT;
LINE_STRUCT = set_defaults_struct(LINE_STRUCT,DEF_LINE_STRUCT);
% PLOT_STRUCT = p.Results.PLOT_STRUCT;
% PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEF_PLOT_STRUCT);
if isempty(LINE_STRUCT.err_edge_color)
    LINE_STRUCT.err_edge_color = 'none';
end
%% ===================================================================== %%
tt = tic;
szx = size(xx);
szy = size(yy);
sze = size(err_bounds);
mx = find(max(szx) == szx);
my = find(max(szx) == szy);
me = find(max(szx) == sze);
if ~all([mx,my,me] == 1)
    dim1 = mx;
    dim2 = setdiff(1:length(szx),mx);
    xx = permute(xx,[dim1,dim2]);
    %--
    dim1 = my;
    dim2 = setdiff(1:length(szy),my);
    yy = permute(yy,[dim1,dim2]);
    %--
    dim1 = me;
    dim2 = setdiff(1:length(sze),me);
    err_bounds = permute(err_bounds,[dim1,dim2]);
end

%##
U = err_bounds(:,2);
L = err_bounds(:,1);
Xcoords = [xx xx(end:-1:1)];
Ycoords = [U L(end:-1:1)];
%## PLOT PATCH & LINE
Pa = patch(Xcoords(:),Ycoords(:),LINE_STRUCT.err_color);
set(Pa,'LineStyle',LINE_STRUCT.err_line_style, ...
    'LineWidth',LINE_STRUCT.err_line_width, ...
    'FaceAlpha',LINE_STRUCT.err_alpha, ...
    'EdgeColor',LINE_STRUCT.err_edge_color);
hold on;
Li = plot(xx,yy, ...
    'LineWidth',LINE_STRUCT.line_width, ...
    'LineStyle',LINE_STRUCT.line_style, ...
    'DisplayName',LINE_STRUCT.line_label);
set(Li,'Color',[LINE_STRUCT.line_color,LINE_STRUCT.line_alpha]);

hold on;
%-- time
fprintf('done. cus_jackknife_sung.m %0.2g s\n',toc(tt));
end

