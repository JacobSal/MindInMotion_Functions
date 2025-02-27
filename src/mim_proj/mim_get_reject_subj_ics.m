function [ic_rej_out] = mim_get_reject_subj_ics(EEG,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%- logo
cat_logo();
%## PATHS
if ~ispc
    pp = strsplit(path,':');
else
    pp = strsplit(path,';');
end
%-
tmp = regexp(pp,'eeglab');
tmp = pp(~cellfun(@isempty,tmp));
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_eeglab = b2{1}(1:b1{1});
    fprintf('EEGLAB path: %s\n',path_eeglab);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('EEGLAB path not found.\n');
    end
end
%## TIME
tt = tic;
%## DEFINE DEFAULTS
DEF_REJ_STRUCT = struct( ...
    'powpow_params',struct('do_calc',true,...
        'upper_freq_lim',100, ...
        'input_data_int',2, ...
        'method_int',2, ...
        'num_iters',300, ...
        'do_plot',true), ...
    'psd_params',struct('fit_freqs',(2:40), ...
        'calc_freqs',[2,100], ...
        'slope_thresh',-0.2, ...
        'spec_perc',80), ...
    'iclabel_params',struct('iclabel_ver','lite', ...
        'class_ints_keep',[1,1,1], ...
        'class_threshs_keep',[0.5,0.75,0.9], ...
        'class_ints_rmv',[2,3], ...
        'class_threshs_rmv',[0.5,0.5], ...
        'class_keep_scores',[8,3,2]),... % NOTE: iclabel classes: 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
    'do_valid_plots',true, ...
    'ica_rv_thresh',0.15, ...
    'plot_save_dir',{''});
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## PARAMETER
addParameter(p,'REJ_STRUCT',DEF_REJ_STRUCT,@(x) validate_struct(x,DEF_REJ_STRUCT));
parse(p,EEG,varargin{:});
%## SET DEFAULTS
REJ_STRUCT = p.Results.REJ_STRUCT;
%-
REJ_STRUCT = set_defaults_struct(REJ_STRUCT,DEF_REJ_STRUCT);
%- other chks
if length(REJ_STRUCT.iclabel_params.class_ints_keep) ~= length(REJ_STRUCT.iclabel_params.class_threshs_keep)
    error('REJ_STRUCT.iclabel_params.class_ints_keep & class_threshs_keep must be same length')
end
if length(REJ_STRUCT.iclabel_params.class_ints_rmv) ~= length(REJ_STRUCT.iclabel_params.class_threshs_rmv)
    error('REJ_STRUCT.iclabel_params.class_ints_rmv & class_threshs_rmv must be same length')
end
%% ===================================================================== %%
%## RECALCULATE ICAACT MATRICES
EEG = eeg_checkset(EEG,'loaddata');
if isempty(EEG.icaact)
    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

% %## (CRITERIA 1) REJECT BASED ON ICLABELS
% params_icl = REJ_STRUCT.iclabel_params;
% try 
%     EEG.etc.ic_classification.ICLabel.classifications;
% catch 
%     EEG = iclabel(EEG,params_icl.iclabel_ver); % NOTE: iclabel classes: 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
% end
% EEG.etc.ic_classification.ICLabel.classifications_default = EEG.etc.ic_classification.ICLabel.classifications;
% iclabel_percs = EEG.etc.ic_classification.ICLabel.classifications;
% %- Count brain components based on ICLabel
% comps_keep = cell(1,length(params_icl.class_threshs_keep));
% comps_remove = cell(2,length(params_icl.class_threshs_rmv));
% %- loop
% for t_i = 1:length(params_icl.class_threshs_keep)
%     tmp = sum(iclabel_percs(:,params_icl.class_ints_keep),2);
%     comps_keep{t_i} = find(tmp>params_icl.class_threshs_keep(t_i));
% end
% %- loop
% for t_i = 1:length(params_icl.class_threshs_rmv)
%     tmp = sum(iclabel_percs(:,2),2); % muscle
%     comps_remove{1,t_i} = find(tmp>params_icl.class_threshs_rmv(t_i));
%     tmp = sum(iclabel_percs(:,3),2); % eye
%     comps_remove{2,t_i} = find(tmp>params_icl.class_threshs_rmv(t_i));
% end
% %-
% % all_ic_crit.ICLabel.brain50 = comps_keep{1};
% % all_ic_crit.ICLabel.brain75 = comps_keep{2};
% % all_ic_crit.ICLabel.brain90 = comps_keep{3};
% %-
% all_ic_crit.ICLabel.brain = comps_keep{1};
% %(01/08/2025) JS, using the 50% cutoff here to be more liberal 
% % all_ic_crit.ICLabel.muscle = comps_remove{1,1};
% % all_ic_crit.ICLabel.muscle_threshold = params_icl.class_threshs_rmv(1);
% % all_ic_crit.ICLabel.eye = comps_remove{2,1};
% % all_ic_crit.ICLabel.eye_threshold = params_icl.class_threshs_rmv(1);
% all_ic_crit.ICLabel.classification = EEG.etc.ic_classification.ICLabel.classifications;
% all_ic_crit.ICLabel.rejection_crit = params_icl;
% %-
% all_ic_crit.Spectra.keep = brain_psd_keep_ics;
% all_ic_crit.Spectra.dump = bad_psd_ics;
% %-
% all_ic_crit.RV.keep = ics_rv_keep;
% all_ic_crit.RV.IC_RV = ic_rv;
% EEG.etc.IC_rej = all_ic_crit;

%% (PSD CRITERIA) ====================================================== %%
params_psd = REJ_STRUCT.psd_params;
%## (Criteria 2) Spectrum graph: Plot PSD together for selected channels
potential_ics = 1:size(EEG.icawinv,2);
tmp_icaact = EEG.icaact(potential_ics, :, :);
%-- compute the spectrum plots: TO DO: need to run it faster
%output unit from spectopo is in db
[spectra_psd,psd_freqs,~,~,~] = spectopo(tmp_icaact, EEG.pnts, EEG.srate, ...
    'percent',params_psd.spec_perc,...
    'freqrange',params_psd.calc_freqs, ...
    'plot','on');

%## CUSTOM PWELCH PSD?
% psd_freq_res = 1;
% psd_freqs = (params_psd.calc_freqs(1):psd_freq_res:params_psd.calc_freqs(2));
% spectra_psd = zeros(size(tmp_icaact,1),length(psd_freqs));
% for ci = 1:size(tmp_icaact,1)
%     padding = repmat(mean(tmp_icaact(ci,:)),[1,ceil(EEG.srate*pi())]);
%     tmp = cat(2,padding,tmp_icaact(ci,:),padding);
%     [Pxx,F] =   pwelch(tmp, ...
%         floor(EEG.srate/pi()),[],length(psd_freqs)*2,EEG.srate);
%     spectra_psd(ci,:) = Pxx'; 
% end
%-- add a linear fit: note that the frequency used here is log(frequency)
%to help characterize the 1/f structure. y axis should be log(power)
lin_fit = zeros(length(potential_ics),2);
spectra_psd_fit = zeros(length(potential_ics),length(psd_freqs));
for i = 1:length(potential_ics)
    lin_fit(i,:) = polyfit(log10(psd_freqs(params_psd.fit_freqs,1)), spectra_psd(i,params_psd.fit_freqs)', 1);
    spectra_psd_fit(i,:) = lin_fit(i,1).*log10(psd_freqs)+lin_fit(i,2);
end
%-- IC rejection criteria: only keep those with log fit < 0
brain_psd_keep_ics = potential_ics(lin_fit(:,1)<params_psd.slope_thresh)';% actually not a lot get discarded
bad_psd_ics = potential_ics(lin_fit(:,1)>=params_psd.slope_thresh)';
num_ics = 1:size(EEG.icawinv,2);

%% (RV CRITERIA) ======================================================= %%
%## (Criteria 4) Scalp topographs and dipole location (if outside the brain or not)
% Residual variance < 0.15
ic_rv = vertcat(EEG.dipfit.model.rv);
ic_posxyz = vertcat(EEG.dipfit.model.posxyz);
ics_rv_keep = find(ic_rv <= REJ_STRUCT.ica_rv_thresh & all(~isnan(ic_posxyz),2));

%% (ICLABEL CRITERIA) ================================================== %%
%## GET ICLABELS
icl_classes = {'brain','muscle','eye','heart','line_noise','channel_noise','other'};
params_icl = REJ_STRUCT.iclabel_params;
try 
    EEG.etc.ic_classification.ICLabel.classifications;
catch 
    EEG = iclabel(EEG,params_icl.iclabel_ver); % NOTE: iclabel classes: 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
end
EEG.etc.ic_classification.ICLabel.classifications_default = EEG.etc.ic_classification.ICLabel.classifications;
iclabel_percs = EEG.etc.ic_classification.ICLabel.classifications;
%-- Count brain components based on ICLabel
comps_keep = cell(1,length(params_icl.class_threshs_keep));
comps_remove = cell(1,length(params_icl.class_threshs_rmv));
%-- keep crit
for t_i = 1:length(params_icl.class_threshs_keep)
    tmp = sum(iclabel_percs(:,params_icl.class_ints_keep(t_i)),2);
    comps_keep{t_i} = find(tmp>params_icl.class_threshs_keep(t_i));
    all_ic_crit.ICLabel.(sprintf('%s%i',icl_classes{params_icl.class_ints_keep(t_i)},params_icl.class_threshs_keep(t_i)*100)) = comps_keep{t_i};
end
%-- removal crit
for t_i = 1:length(params_icl.class_ints_rmv)
    tmp = sum(iclabel_percs(:,params_icl.class_ints_rmv(t_i)),2);
    comps_remove{t_i} = find(tmp>params_icl.class_threshs_rmv(t_i));
    % all_ic_crit.ICLabel.(sprintf('%s%i',icl_classes{params_icl.class_ints_rmv(t_i)},params_icl.class_threshs_rmv(t_i)*100)) = comps_remove{t_i};
    all_ic_crit.ICLabel.(sprintf('%s',icl_classes{params_icl.class_ints_rmv(t_i)})) = comps_remove{t_i};
end

%% (CRITERIA WEIGHTING) ================================================ %%
%## GATHER ICLABEL DATA & OTHER CRITERIA
all_ic_crit.ICLabel.brain = comps_keep{1};
%(01/08/2025) JS, using the 50% cutoff for brain to be more liberal
all_ic_crit.ICLabel.classification = EEG.etc.ic_classification.ICLabel.classifications;
all_ic_crit.ICLabel.rejection_crit = params_icl;
%--
all_ic_crit.Spectra.keep = brain_psd_keep_ics;
all_ic_crit.Spectra.dump = bad_psd_ics;
%--
all_ic_crit.RV.keep = ics_rv_keep;
all_ic_crit.RV.IC_RV = ic_rv;
EEG.etc.IC_rej = all_ic_crit;

%## ASSIGN WEIGHTS TO CRITERIA
all_eye_ics = zeros(size(EEG.icawinv,2),1);
all_muscle_ics = zeros(size(EEG.icawinv,2),1);
all_brain_ics = zeros(size(EEG.icawinv,2),1);
%-- If muscle higher score remove
all_muscle_ics(all_ic_crit.ICLabel.muscle) = all_muscle_ics(all_ic_crit.ICLabel.muscle)+2;
all_muscle_ics(all_ic_crit.Spectra.dump) = all_muscle_ics(all_ic_crit.Spectra.dump)+1;
%     IC_all_muscle(All_IC_criteria.Projection.EMG) = IC_all_muscle(All_IC_criteria.Projection.EMG)+2;

%-- If eye higher score remove
all_eye_ics(all_ic_crit.ICLabel.eye) = all_eye_ics(all_ic_crit.ICLabel.eye)+2;
all_eye_ics(all_ic_crit.Spectra.dump) = all_eye_ics(all_ic_crit.Spectra.dump)+1;

%-- If brain higher score keep
all_brain_ics(all_ic_crit.ICLabel.brain50) = all_brain_ics(all_ic_crit.ICLabel.brain50)+2;
all_brain_ics(all_ic_crit.ICLabel.brain75) = all_brain_ics(all_ic_crit.ICLabel.brain75)+2;    
all_brain_ics(all_ic_crit.Spectra.keep) = all_brain_ics(all_ic_crit.Spectra.keep)+1;
%     IC_all_brain(All_IC_criteria.Projection.EMG) = IC_all_brain(All_IC_criteria.Projection.EMG)-3;
all_brain_ics(all_ic_crit.RV.keep) = all_brain_ics(all_ic_crit.RV.keep)+5;
all_brain_ics(size(EEG.icawinv,2)-5:size(EEG.icawinv,2)) = all_brain_ics(size(EEG.icawinv,2)-5:size(EEG.icawinv,2))-1;
%(01/08/2025) JS, if they are really low on the ICA "list" they probably didn't seperate
%well and are therefore less likely to be brain. Using last 5 as crit
%% (POWPOWCAT CRITERIA) ================================================ %%
%## (Criteria 5)  PowPowCat
% PowPow Cat Cross-Frequency Power-Power Coupling Analysis: A Useful Cross-Frequency Measure to Classify ICA-Decomposed EEG
% It takes a long time to run. should pick only the ones that are
% classified as 'brain'
bad_powpow_ics = [];
params_pp = REJ_STRUCT.powpow_params;
if params_pp.do_calc
%     IC_powpow = find(Output_ICRejection.IC_all_brain >= 8);
    fprintf('PowPowCAT parameters:\n upperFreqLimit= %i Hz\n inputDataType = ICs\n methodType= Spearman''s correlation (non-parametric)\n numIterations = %i\n', ...
        params_pp.upper_freq_lim,params_pp.num_iters);
    %run PowPowCAT
    powpow_ics = find(all_brain_ics >= params_icl.class_keep_scores(1)-1);
%     IC_powpow = find(IC_all_brain >= 5);
    if ~isempty(powpow_ics)
        tmp = EEG;
        tmp.icaact = EEG.icaact(powpow_ics,:);
        fprintf('\n');
        %## CALC POWPOW
        tmp = calc_PowPowCAT(tmp, ...
            params_pp.upper_freq_lim, ...
            params_pp.input_data_int, ...
            params_pp.method_int, ...
            params_pp.num_iters);
        fprintf('\n');
        tmp.setname = 'powpowcat';        
        %## USE POWPOW CRIT TO REJECT ICS
        bad_ics_out = PowPowCat_ICrej(tmp,params_pp.do_plot,REJ_STRUCT.plot_save_dir);
        %## POP_PROPS_EXTENDED
        if params_pp.do_plot
            tmp_save_dir = [REJ_STRUCT.plot_save_dir filesep 'powpowcat'];
            if ~exist(tmp_save_dir,'dir')
                mkdir(tmp_save_dir)
            end
            
            %## (PLOT 1) POWPOWCAT ALL IC PLOT
            Plot35IcPushbutton_powpow(tmp,length(powpow_ics),powpow_ics)
            fig_i = get(groot,'CurrentFigure');
            exportgraphics(fig_i,[tmp_save_dir filesep 'powpowcat.jpg']);

            %## (PLOT 2) EXTENDED COMPONENT PROPERTIES
            for ii = 1:length(bad_ics_out)
                fprintf('\nPlotting extended component properties of bad ICs');
                % plot extended comp properties
                pop_prop_extended(tmp,0,bad_ics_out(ii),nan(),{'freqrange',[2,80]},{},1,'')
                fig_i = get(groot,'CurrentFigure');
                exportgraphics(fig_i,[tmp_save_dir filesep sprintf('%s_properties_%i.jpg',tmp.subject,powpow_ics(bad_ics_out(ii)))], ...
                    'Resolution',300);
                close(fig_i);
            end           
        end
        bad_powpow_ics = powpow_ics(bad_ics_out);
    end
end
%% (SAVE OUTPUT) ======================================================= %%
ic_rej_out.All_IC_criteria = all_ic_crit;
ic_rej_out.IC_all_brain = all_brain_ics;
ic_rej_out.IC_all_muscle = all_muscle_ics;
ic_rej_out.IC_all_eye = all_eye_ics;
ic_rej_out.IC_powpow_rej = bad_powpow_ics;
ic_rej_out.Cleaning_Params = EEG.etc.Params;
%% (VALIDATION PLOTS) ================================================== %%
%## (PLOTS) INSPECTION
if REJ_STRUCT.do_valid_plots
    %## (PLOT 2) STEMPLOT OF PSD CRITERIA SCORES
    figure();
    hold on;
    stem(num_ics,lin_fit(:,1));
    if any(lin_fit(:,1)> params_psd.slope_thresh)
        stem(num_ics(lin_fit(:,1)> params_psd.slope_thresh),lin_fit(lin_fit(:,1)>params_psd.slope_thresh,1),'r');
    end
    hold off;
    ylabel('slope')
    fig_i = get(groot,'CurrentFigure');
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep sprintf('spectral_stem.jpg')]);
    close(fig_i);

    %## (PLOT 3) GOOD ICS PSD PLOTS
    % Log-log plot - retained
    cmap_in = linspecer(length(brain_psd_keep_ics)); % my personal color scheme 
    fig_i = figure('color','w');
    subplot(1,2,1)
    hold on;
    for i0 = 1:length(brain_psd_keep_ics)
        % plot the fitted line
        p = brain_psd_keep_ics(i0);
        semilogx(psd_freqs,lin_fit(p,1).*log10(psd_freqs)+lin_fit(p,2),'--','color',cmap_in(i0,:),'linewidth',1.2);hold on;
        semilogx(psd_freqs,spectra_psd(p,:),'color',cmap_in(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title('KEEP: IC Activity Power Spectrum','units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(brain_psd_keep_ics)),'Location','eastoutside');
    legend box off
    subplot(1,2,2)
    hold on;
    for i0 = 1:length(brain_psd_keep_ics)
        % plot the fitted line
        p = brain_psd_keep_ics(i0);
        plot(psd_freqs,spectra_psd(p,:),'color',cmap_in(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title('KEEP: IC Activity Power Spectrum','units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(brain_psd_keep_ics)),'Location','eastoutside');
    legend box off
    % fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [0 0 1080 1080];
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep 'keep_potential_brain_components_spectral.jpg']);
    close(fig_i);

    %## (PLOT 4) BAD ICS PSD PLOTS
    cmap_in = linspecer(length(bad_psd_ics));%my personal color scheme 
    figure('color','w');
    hold on;
    for i0 = 1:length(bad_psd_ics)
        % plot the fitted line
        p = bad_psd_ics(i0);
        semilogx(psd_freqs,lin_fit(p,1).*log10(psd_freqs)+lin_fit(p,2),'--','color',cmap_in(i0,:),'linewidth',1.2);hold on;
        semilogx(psd_freqs,spectra_psd(p,:),'color',cmap_in(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title('BAD: IC Activity Power Spectrum','units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(bad_psd_ics)),'Location','eastoutside');
    legend box off
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [0 0 1080 720];
    %     saveas(gcf,fullfile(save_IC_Rejection_folder,subjStr,'Figures',['potential brain components_spectral_',subDirNum,'.fig']))
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep 'reject_potential_brain_components_spectral.jpg'])
    close(fig_i);

    %## (PLOT 5) STEMPLOT OF ICLABEL SCORES
    fig_i = figure('color','white');
    hold on;
    subplot(3,1,1);
    stem(num_ics',all_muscle_ics);
    title(['muscle: ',num2str(sum(all_muscle_ics == params_icl.class_keep_scores(2)))]);
    ylabel('Dump: score');
    subplot(3,1,2);
    stem(num_ics',all_eye_ics);
    title(['eye: ',num2str(sum(all_eye_ics == params_icl.class_keep_scores(3)))]);
    ylabel('Dump: score');
    hold on;
    ss = subplot(3,1,3);    
    stem(num_ics',all_brain_ics);
    title(['brain: ',num2str(sum(all_brain_ics == params_icl.class_keep_scores(1)))]);
    ylabel('Keep: score');
    stem(num_ics(all_muscle_ics>=params_icl.class_keep_scores(2)-1)', ...
        all_brain_ics(all_muscle_ics>=params_icl.class_keep_scores(2)-1), ...
        'r', ...
        'DisplayName','flagged muscle');
    xlabel('IC number');
    legend(ss);
    % legend('','flag muscle');
    % fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [0 0 1080 720];
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep 'IC_score.jpg'])
    close(fig_i);

    %## (PLOT 6) LOG-LOG PLOT 
    % Log-log plot - retained
    powpow_ics = find(all_brain_ics >= params_icl.class_keep_scores(1)-1);
    cmap_in = linspecer(length(powpow_ics));
    fig_i = figure('color','w');
    subplot(1,2,1)
    for i0 = 1:length(powpow_ics)
        % plot the fitted line
        p = powpow_ics(i0);
        semilogx(psd_freqs,lin_fit(p,1).*log10(psd_freqs)+lin_fit(p,2),'--','color',cmap_in(i0,:),'linewidth',1.2);hold on;
        semilogx(psd_freqs,spectra_psd(p,:),'color',cmap_in(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    title('KEEP: IC Activity PSD','units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('PSD (dB)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(reshape(repmat(cellstr(num2str(powpow_ics)),2)',[],1),'Location','eastoutside');
    legend box off
    subplot(1,2,2)
    for i0 = 1:length(powpow_ics)
        % plot the fitted line
        p = powpow_ics(i0);
        plot(psd_freqs,spectra_psd(p,:),'color',cmap_in(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    title('KEEP: IC Activity PSD','units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('PSD (dB)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(powpow_ics)),'Location','eastoutside');
    legend box off
    fig_i.Position = [0 0 1080 720];
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep 'KEEP_IC_PSD.jpg']);
    close(fig_i);

    %## (PLOT 7) BEMOBIL PIPELINE PLOT
    summed_scores_to_keep = sum(iclabel_percs(:,1),2);
    titles_FEM = cell(size(EEG.icawinv,2),1);
    for i_title = 1:size(EEG.icawinv,2)
        titles_FEM{i_title} = [num2str(i_title) ': ' num2str(round(summed_scores_to_keep(i_title),3))];
    end
    bemobil_plot_patterns_CL(EEG.icawinv,EEG.chanlocs,'chan_to_plot',find(all_brain_ics >= 8)','weights',summed_scores_to_keep,...
        'fixscale',0,'minweight',0.75,'titles',titles_FEM);
    fig_i = get(groot,'CurrentFigure');
    exportgraphics(fig_i,[REJ_STRUCT.plot_save_dir filesep 'potential_brain_components_allcritera_customElectrode.jpg']);
    close(fig_i);
end
fprintf('done. mim_get_reject_subj_ics.m: %0.2g s\n',toc(tt));
end

