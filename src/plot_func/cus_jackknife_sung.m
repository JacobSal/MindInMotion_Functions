function [ax,Pa,Li] = cus_jackknife_sung(ax,xx,yy,err_bounds,varargin)
% JACKKNIFE plots jackknife errorbars around a given curve
%
%     Plots a gray Jackknife around the line displayed in black very useful for
%     funky nature style error bars which are shaded.
%
%         Pa is a patch object for more help on patch objects see below
%         Li is a line object, more help on line object is available in MATLAB
%
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
PRINT_CHKS = false;
%--
DEF_LINE_STRUCT = struct('do_line_avg',false, ...
    'line_props',{{ ...
        'LineWidth',2, ...
        'LineStyle','-', ...
        'DisplayName','line', ...
        }}, ...
    'line_color',[0,0,0,1], ...
    'line_avg_fcn',@(x) mean(x,1), ...
    'do_err_shading',true, ...
    'err_props',{{ ...
        'LineStyle',':', ...
        'LineWidth',3, ...
        'FaceAlpha',0.6, ...
        'EdgeColor','none', ...
        'FaceColor',[0.5,0.5,0.5]}}, ...
    'err_bnd_vec',[]);

err_bnd_struct = struct( ...
    'LineStyle',':', ...
    'LineWidth',3, ...
    'FaceAlpha',0.6, ...
    'EdgeColor','none', ...
    'FaceColor',[0.5,0.5,0.5]);
psd_line_struct = struct('LineWidth',2, ...
    'LineStyle','-', ...
    'DisplayName','line', ...
    'Color',[0.5,0.5,0.5,0.65] ...
    );
DEF_LINE_STRUCT = struct('do_line_avg',false, ...
    'line_props',psd_line_struct, ...
    'line_avg_fcn',@(x) mean(x,1), ...
    'do_err_shading',true, ...
    'err_bnd_props',err_bnd_struct, ...
    'err_bnd_vec',[]);

%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'ax',@isgraphics);
addRequired(p,'xx',@isnumeric);
addRequired(p,'yy',@isnumeric);
addRequired(p,'err_bounds',@(x) isnumeric(x) && any(size(x) == 2));
%-- PARAMETER
addParameter(p,'LINE_STRUCT',DEF_LINE_STRUCT,@(x) validate_struct(x,DEF_LINE_STRUCT,PRINT_CHKS));
parse(p,ax,xx,yy,err_bounds,varargin{:});
%-- SET DEFAULTS
LINE_STRUCT = p.Results.LINE_STRUCT;
LINE_STRUCT = set_defaults_struct(LINE_STRUCT,DEF_LINE_STRUCT,PRINT_CHKS);

%## HANDLE COLOR MAPPING
if isempty(LINE_STRUCT.line_props.Color) && ~isempty(LINE_STRUCT.line_props.EdgeColor)
    LINE_STRUCT.line_props.Color = LINE_STRUCT.line_props.EdgeColor;
end

line_props = struct2args(LINE_STRUCT.line_props);
err_bnd_props = struct2args(LINE_STRUCT.err_bnd_props);
%% ===================================================================== %%
tt = tic;
%##
U = err_bounds(:,2);
L = err_bounds(:,1);
Xcoords = [xx xx(end:-1:1)];
Ycoords = [U L(end:-1:1)];
%## PLOT PATCH & LINE
Pa = patch(Xcoords(:),Ycoords(:),'w', ...
    err_bnd_props{:});
hold on;
Li = plot(xx,yy, ...
    line_props{:});
hold on;
%-- time
fprintf('done. cus_jackknife_sung.m %0.2g s\n',toc(tt));
end

