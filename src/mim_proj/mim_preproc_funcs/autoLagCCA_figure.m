% autoCCA_figure
%% visualize spectra results - What is it visualizing? 
function autoLagCCA_figure(EEG)
    visCCA_PSD = 1;
    if visCCA_PSD
        EEG_preCleaning = EEGandEMG;
        EEG_postCleaning = EEG_CCA;
        xChanPlot = 1:EEG_preCleaning.nbchan;
        % timeWindow = find(EEG_preCleaning.times/1000 >= 0 & EEG_preCleaning.times/1000 <= 200); %arbitrarily calculating PSD no smaller window so we don't have to wait forever
        timeWindow = find(EEG_preCleaning.times/1000 >= 0 & EEG_preCleaning.times/1000 <= inf);
    %         [spectra_in,FREQ,compeegspecdB_in,resvar_in,specstd_in] = spectopo(EEG_preCleaning.data(xChanPlot,timeWindow), length(timeWindow), EEG_preCleaning.srate, 'freqfac',2, 'plot','off');
    %         [spectra_out,FREQ,compeegspecdB_out,resvar_out,specstd_out] = spectopo(EEG_postCleaning.data(xChanPlot,timeWindow), length(timeWindow), EEG_postCleaning.srate, 'freqfac',2, 'plot','off');
        tempDiffMat =EEG_postCleaning.data(xChanPlot,timeWindow)-EEG_preCleaning.data(xChanPlot,timeWindow);
        if any(tempDiffMat(:)~=0) %make sure there is some non-zero difference between the datasets somewhere in any channel or time point
            [spectra_diff,FREQ] =  spectopo( tempDiffMat, length(timeWindow), EEG.srate, 'freqfac',1, 'plot','off');
        else
            spectra_diff = nan(size(spectra_in)); disp('No difference between the clean and dirty data sets!')
        end
        % rnk = rank(double(tempEEG.data)); % Calculate rank of output data (component removal reduces rank)
        % fprintf(['CCA EEG rank out: ' int2str(rnk) '\n'])

        figure; % Plot spectra of mean channel input and output
        set(gcf,'Position',[500 200 500 400]);
        set(gcf,'Color','w');
        title('Power spectra');
        hold on;
    %         plot(FREQ, median(spectra_in,1), 'r','LineWidth', .5); plot(FREQ, prctile(spectra_in,25), 'r--','LineWidth', .5); plot(FREQ, prctile(spectra_in,75), 'r--','LineWidth', .5);
    %         set(gca,'FontSize',12);
    %         plot(FREQ, median(spectra_out,1), 'b','LineWidth', .5); plot(FREQ, prctile(spectra_out,25), 'b--','LineWidth', .5); plot(FREQ, prctile(spectra_out,75), 'b--','LineWidth', .5);
        plot(FREQ, median(spectra_diff,1), 'k','LineWidth',.5);
        legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Diff Q2','Location', 'eastoutside', 'Orientation', 'vertical');
        legend boxoff;
        xlabel('Frequency (Hz)');
        ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
        hold off;
        xlim([.5 70]);  set(gca,'Box','off');
        YLIM = get(gca,'YLim');
        set(gca,'XTick', [0.5,4,8,13,30,60]);
        set(gca, 'XTickLabel', {'','4','8','13','30','60'});
%         saveas(gcf,strcat(fullfile(MiM_config.MIMDataFolder_CL,subjStr,'EEG','Figures'),'\PSD post-CCA','.jpg'));
        %end visualizing PSD
    end
end