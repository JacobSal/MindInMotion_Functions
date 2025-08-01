function y=gety(ax,varargin)
% Helper function that Returns the largest single value of ydata in a given
% graphic handle. It returns the given value if it is not a graphics
% handle. 
% Returns the largest single value of ydata in a given graphic handle
% h= figure,axes,line. Note that y=h if h is not a graphics
%## TIME
tic
%## DEFINE DEFAULTS
%- this preference depends on graphics type being plotted
% graphic_types = {'line','hggroup','patch'}; %,'scatter'
graphic_types = {'line','hggroup','patch','scatter'}; 
%## PARSER
p = inputParser;
%- REQUIRED
addRequired(p,'ax')
addOptional(p,'graphic_types',graphic_types,@iscell)
parse(p,ax,varargin{:});
%% ===================================================================== %%
if ~isempty(ax)
    switch get(ax,'type')
        case graphic_types %{'line','hggroup','patch','scatter'},
            y = max(get(ax,'ydata'));
            return;
        otherwise
            ys = [];
            hs = get(ax,'children');
            for n = 1:length(hs)
                ys = [ys,gety(hs(n))];
            end
            y = max(ys(:));
    end
    if isempty(y)
        try 
            tmp = ax.Units;
            if ~strcmp(tmp,'data')
                ax.Units = 'data';
                y = ax.Position(2);
                ax.Units = tmp;
            else
                y = ax.Position(2);
            end
        catch e
            fprintf('ax is not a position element.\n%s\n\n',getReport(e));
            y = 0;
        end
    end
else
    y = 0;
end
toc
end

