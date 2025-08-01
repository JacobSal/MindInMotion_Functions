function [struct_out] = set_defaults_struct(x,DEFAULT_STRUCT,varargin)
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
struct_out = x;
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
s_flags = zeros(length(fs2),1);
vals1 = struct2cell(x);
%-- check field value's class type
pchk_fun(sprintf('Setting structure defaults...\n'),print_chks)
for f = 1:length(fs2)
    if isstruct(DEFAULT_STRUCT.(fs2{f}))
        s_flags(f) = true;
    end
end
%-- 
for f = 1:length(fs2)
    ind1 = strcmp(fs2{f},fs1);
    if any(ind1)
        if isempty(vals1{ind1}) && ~ischar(vals1{ind1})
            pchk_fun(sprintf('Setting struct.%s to default value\n',fs2{f}),print_chks)
        
            struct_out.(fs1{ind1}) = DEFAULT_STRUCT.(fs2{f});
        elseif isstruct(DEFAULT_STRUCT.(fs2{f}))
            pchk_fun(sprintf('Recursive iteration on struct.%s...\n',fs2{f}),print_chks)
        
            tmp = set_defaults_struct(x.(fs2{f}),DEFAULT_STRUCT.(fs2{f}),print_chks);
            struct_out.(fs1{ind1}) = tmp;
        end
    else
        pchk_fun(sprintf('Setting struct.%s to default value\n',fs2{f}),print_chks)
        struct_out.(fs2{f}) = DEFAULT_STRUCT.(fs2{f});
    end
end
end

%##
function pchk_fun(str_in,print_chks)
    if print_chks
        fprintf(str_in);
    end
end

