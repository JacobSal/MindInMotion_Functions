function [] = par_save(SAVEVAR,fpath,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           SAVEVAR, variable to save.
%           fPath, CHAR
%               path to the folder or file
%       OPTIONAL:
%           fName, CHAR | EMPTY
%               file name & extension (e.g., 'INEEG.mat')
%       PARAMETER:
%   OUT:
%       NONE
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/06/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tt = tic;
%## DEFINE DEFAULTS
%- fPath
errorMsg = 'Value ''fpath'' must be CHAR. The path must exist.'; 
fp_validFcn = @(x) assert(ischar(x),errorMsg);
%-
fname = [];
errorMsg = 'Value ''fname'' must be CHAR || EMPTY.';
fn_validFcn = @(x) assert(ischar(x) || isempty(x),errorMsg);
%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'SAVEVAR')
addRequired(p,'fpath',fp_validFcn)
addOptional(p,'fname',fname,fn_validFcn)
%## PARAMETER
parse(p, SAVEVAR, fpath, varargin{:});
fname = p.Results.fname;
%% ===================================================================== %%
% set save path
%- if filename is included in path
if isempty(fname)
    spath = strsplit(fpath,filesep);
    jpath = join(spath(1:end-1)',filesep);
    fname = spath{end};
    fpath = jpath{1};  
end
%-
save_fpath = [fpath filesep fname];
if ~exist(fpath,'dir')
    mkdir(fpath)
end
%- save
s = whos('SAVEVAR');
fprintf('\n%s is %0.2g Gigabytes\n',fname,s.bytes/1e9);
if s.bytes >= 2e9
    fprintf('Saving %s using ''v7.3'' to\n%s\n',fname,fpath);
    save(save_fpath,'SAVEVAR','-v7.3');
else
    fprintf('Saving %s using ''v6'' to\n%s\n',fname,fpath);
    save(save_fpath,'SAVEVAR','-v6');
end
fprintf('done. par_save.m: %0.2f s\n\n',toc(tt));
end

