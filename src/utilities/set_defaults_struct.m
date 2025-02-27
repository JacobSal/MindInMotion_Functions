function [struct_out] = set_defaults_struct(x,DEFAULT_STRUCT)
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
struct_out = x;
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
s_flags = zeros(length(fs2),1);
vals1 = struct2cell(x);
%- check field value's class type
fprintf('Setting structure defaults...\n');
for f = 1:length(fs2)
    if isstruct(DEFAULT_STRUCT.(fs2{f}))
        s_flags(f) = true;
    end
end
for f = 1:length(fs2)
    ind1 = strcmp(fs2{f},fs1);
    if any(ind1)
        if isempty(vals1{ind1}) && ~ischar(vals1{ind1})
            fprintf(2,'Setting struct.%s to default value\n',fs2{f});
            struct_out.(fs1{ind1}) = DEFAULT_STRUCT.(fs2{f});
        elseif isstruct(DEFAULT_STRUCT.(fs2{f}))
            fprintf('Recursive iteration on struct.%s...\n',fs2{f});
            tmp = set_defaults_struct(x.(fs2{f}),DEFAULT_STRUCT.(fs2{f}));
            struct_out.(fs1{ind1}) = tmp;
        end
    else
        fprintf(2,'Setting struct.%s to default value\n',fs2{f});
        struct_out.(fs2{f}) = DEFAULT_STRUCT.(fs2{f});
    end
end
end
