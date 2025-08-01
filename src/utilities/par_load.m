function [OUTVAR] = par_load(fpath,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           fpath, CHAR
%               path to the folder where your file is held
%           fname, CHAR
%               file name & extension (e.g., 'INEEG.mat')
%       OPTIONAL:
%           path_ext, CHAR (default: [])
%               extension for operating system conversion see.
%               convertPath2Drive.m & convertPath2UNIX.m
%       PARAMETER:
%   OUT:
%       OUTVAR, VAR
%           
%   IMPORTANT:
%## TIME
tt = tic();
%## DEFINE DEFAULTS
fname = [];
errorMsg = 'Value ''fname'' must be CHAR || EMPTY.';
fn_validFcn = @(x) assert(ischar(x) || isempty(x),errorMsg);
struct_fn = [];

%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'fpath',@ischar)
%## OPTIONAL
addOptional(p,'fname',fname,fn_validFcn)
addOptional(p,'struct_fn',struct_fn,@(x) ischar(x) || isempty(x))
%## PARSE
parse(p, fpath, varargin{:});
%## SET DEFAULTS
fname = p.Results.fname;
struct_fn = p.Results.struct_fn;
%% ===================================================================== %%
%## SET LOAD FPATH
%- if filename is included in path
if isempty(fname)
    spath = strsplit(fpath,filesep);
    jpath = join(spath(1:end-1)',filesep);
    fname = spath{end};
    fpath = jpath{1};  
end

%## LOAD
load_fpath = [fpath filesep fname];
%-- 
s = dir(load_fpath);
if isempty(s)
    error('File %s does not exist.',fname);
end
fprintf('\n%s is %0.2g Gigabytes\n',fname,s.bytes/1e9);
OUTVAR = load(load_fpath);
OUTVAR = OUTVAR.SAVEVAR;
if ~isempty(struct_fn)
    OUTVAR = OUTVAR.(struct_fn);
end
%-- time
fprintf('done. par_load duration: %0.2f s\n',toc(tt))
end

