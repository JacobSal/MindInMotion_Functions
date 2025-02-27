function [table_out] = util_resolve_table(table_in,varargin)
%EEGLAB_RESOLVE_STRUCT Summary of this function goes here
%
%   IN:
%       table_in, CELL ARRAY of STRUCTS
%           cell array of structures that you want to concatentate, but
%           they have mismatched fields. The algorithm will add missing
%           fields to indices lacking them.
%       Optional; struct_ids, CELL ARRAY of CHARS
%           cell array of characters that can be used for verbose debugging
%           purposes. (e.g., cell array of subject names).
%   OUT:
%       struct_out, STRUCT ARRAY
%
%   Version History --> See details at the end of the script.
%   Previous Version: n/a
%   Summary:  
%
table_ids = {};
struct_ids_vfcn = @(x) iscell(x) && length(x) == length(table_in);
cat_logo();
%## TIME
t = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'table_in',@iscell);
%## OPTIONAL
addOptional(p,'struct_ids',table_ids,struct_ids_vfcn)
%## PARSE
parse(p, table_in, varargin{:});
%## SET DEFAULTS
table_ids = p.Results.struct_ids;
%% ===================================================================== %%
table_in = table_in(~cellfun(@isempty,table_in));
% table_fields = cell(1,length(table_in));
def_table_log = struct('ind',0, ...
    'var_names',{''}, ...
    'types',{''});
table_log = repmat(def_table_log,[length(table_in),1]);
for i = 1:length(table_in)
    table_log(i).var_names = table_in{i}.Properties.VariableNames; %fields(table_in{i})';
    table_log(i).types = varfun(@class,table_in{i},'OutputFormat','cell');
    table_log(i).ind = i;
end
table_fields = unique([table_log.var_names]);
iter_fields = table_fields;
for i = 1:length(table_in)
    tmp_table = table_in{i};
    tmp_fields = tmp_table.Properties.VariableNames; %fields(tmp_struct);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,tmp_fields)),iter_fields,'UniformOutput',false); 
    out = [out{:}];
    var_add = iter_fields(~out);
    if any(~out)
        for j = 1:length(var_add)
            tmp_table.(var_add{j}) = zeros(size(tmp_table,1),1);
            if isempty(table_ids)
                fprintf('Index %i) Adding fields %s\n',i,var_add{j});
            else
                fprintf('%s) Adding fields %s\n',table_ids{i},var_add{j});
            end
        end
    end 
%     table_in{subj_i} = EEG;
    % table_in{i} = orderfields(tmp_table);
    table_in{i} = tmp_table;
end
%- CONCATENATE table_in
table_out = cat(1,table_in{:}); %cellfun(@(x) [[]; x], table_in);
fprintf('util_resolve_struct processing time: %0.2f\n\n',toc(t));
end

