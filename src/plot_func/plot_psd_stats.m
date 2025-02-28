% -----------------
% highlight regions
% Function extracted from EEGLab
% -----------------
function [axsignif,tmph] = plot_psd_stats(ax, times, regions, highlightmode, myxlabel)
% if isempty(patchcolor)
%     patchcolor = [0.75 0.75 0.75];
% end
tmph = [];
patchcolor = [0.75 0.75 0.75];
color2 = [0 0 0];
orig_yl = ax.YLim;
y1  = ax.YLim; 
% yl(1) = yl(1)-max(abs(yl));
% yl(2) = yl(2)+max(abs(yl));
if y1(1) > 0
    y1(1) = y1(1)-(5*abs(y1(1)));
else
    y1(1) = y1(1)*5;
end
if y1(2) > 0
    y1(2) = y1(2)*5;
else
    y1(2) = y1(2)+(abs(y(1))*5);
end
if ~strcmpi(highlightmode, 'background')
    pos = get(ax, 'position');
    set(gca, 'xtick', []);
    axsignif = axes('position', [pos(1) pos(2)-pos(4)*0.01 pos(3) pos(4)*0.01 ]);
%     plot(times, times, 'w');
    set(axsignif, 'ytick', []);
    yl2 = ax.YLim;
    yl2(1) = yl2(1)-max(abs(yl2));
    yl2(2) = yl2(2)+max(abs(yl2));
    xlim([times(1) times(end)]);
    xlabel(myxlabel);
    
else
    axsignif = [];
    xlabel(myxlabel);
end

if ~isempty(regions)
%     axes(ax);
    in_a_region = 0;
    for index=1:length(regions)
        if regions(index) && ~in_a_region
            tmpreg(1) = times(index);
            in_a_region = 1;
        end
        if (~regions(index) || index == length(regions)) && in_a_region
            tmpreg(2) = times(min(length(times), index));
            in_a_region = 0;
            if strcmpi(highlightmode, 'background') %indicate significance in the background
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [y1(1) y1(1) y1(2) y1(2)], patchcolor); hold on;
                set(tmph, 'edgecolor', 'none','facealpha',0.5,'edgealpha',0.2);
            else
                oldax = ax;
                axes(axsignif);
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl2(1) yl2(1) yl2(2) yl2(2)], color2); hold on;
                set(tmph, 'edgecolor', color2);
                axes(oldax);
            end
        end
    end
    ylim(orig_yl);
end
  