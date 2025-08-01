function [b] = validate_struct(x,DEFAULT_STRUCT,varargin)
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
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'x')
addRequired(p,'DEFAULT_STRUCT');
%## OPTIONAL
addOptional(p,'print_chks',true,@islogical);
%##
parse(p,x,DEFAULT_STRUCT,varargin{:});
print_chks = p.Results.print_chks;
%% ===================================================================== %%
b = false;
struct_name = inputname(2);
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
vals1 = struct2cell(x);
vals2 = struct2cell(DEFAULT_STRUCT);
%-- check field value's class type
pchk_fun(sprintf('Running structure checks...\n'),print_chks)
%-- check field names
chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
if ~all(chk)
    pchk_fun(sprintf('Remove field(s): %s\n',strjoin(fs1(~chk),',')),print_chks)
    pchk_fun(sprintf('Fields for struct do not match for %s\n',struct_name),print_chks)
    return
end
%- check field value's class type
for f = 1:length(fs2)
    ind = strcmp(fs2{f},fs1);
    if any(ind)
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            pchk_fun(sprintf('\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind})),print_chks)
            return
        end
    else
        pchk_fun(sprintf('Struct.%s is not present in input structure. Setting to default.\n',fs2{f}),print_chks)
    end
end
b = true;
end

%##
function pchk_fun(str_in,print_chks)
    if print_chks
        fprintf(str_in);
    end
end

