function bemobil_plot_EEG_visartifact(EEG_new,EEG_old)
% plot EEG artifact comparison 
% Chang Liu - adapt from bemobil pipeline
% TO DO: make this script run faster...
    plotfigure = figure('color','w');
    set(plotfigure, 'Position', get(0,'screensize'))
    ax1 = subplot(231);
    ax2 = subplot(232);
    ax3 = subplot(233);
    ax4 = subplot(234);
    ax5 = subplot(235);
    ax6 = subplot(236);

    starttime = EEG_old.times(end)/7*1;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax1,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 1 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax1);
    close(fighandle)

    starttime = EEG_old.times(end)/7*2;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax2,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 2 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax2);
    close(fighandle)

    starttime = EEG_old.times(end)/7*3;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax3,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 3 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax3);
    close(fighandle)

    starttime = EEG_old.times(end)/7*4;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax4,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 4 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax4);
    close(fighandle)

    starttime = EEG_old.times(end)/7*5;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax5,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 5 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax5);
    close(fighandle)

    starttime = EEG_old.times(end)/7*6;
    vis_artifacts(EEG_new,EEG_old,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000),...
        'ChannelSubset',find(strcmpi('EEG',{EEG_old.chanlocs.type}))); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax6,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 6 of ' num2str(round(EEG_old.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax6);
    close(fighandle)
end
