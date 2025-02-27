function [outputArg1,outputArg2] = proc_kin_psd_meas(EEG,eeg_icaspec,varargin)
%EXTRACT_GAIT_EVENTS Summary of this function goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 01/10/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%##
tic
%% DEFINE DEFAULTS
if ~which('quaternRotate')
    error('You need to add the Gait-Tracking-With-x-IMU-master free library to get quaternion math');
end
HPfilter_passBand   = 0.2; %Hi Colton. Ryan here. I decided to with 0.2 Hz instead of 0.25. Let me know how it looks
LPfilter_passBand   = 20; %Hz
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'eeg_icaspec',@eeg_icaspec);
%## OPTIONAL
%## PARAMETER
parse(p, EEG, eeg_icaspec, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%% ===================================================================== %%
nolog_eeg_psd = 10.^(eeg_psd/10);
freqs = tmp.freqs;
trialinfo = tmp.trialinfo;

%## RUN FOOOF
fprintf('Running FOOOF...\n');
settings = struct('peak_width_limits',[1,8],...
    'min_peak_height',0.05,...
    'max_n_peaks',5);
f_range = [3, 40];
tmp_f = find(freqs > f_range(1) & freqs < f_range(2));
tmp_f = sort([tmp_f; min(tmp_f)-1;max(tmp_f)+1]);
return_model = true;
tmp = zeros(length(tmp_f),size(eeg_psd,2),size(eeg_psd,3));
for j = 1:size(nolog_eeg_psd,3)
    fr = fooof(freqs,mean(squeeze(nolog_eeg_psd(:,:,j)),2),f_range,settings,return_model);
    for i = 1:size(nolog_eeg_psd,2)
        % fr = fooof(freqs,mean(squeeze(nolog_eeg_psd(:,:,j)),2),f_range,settings,return_model);
        % fooof_diff = 10*(fr.power_spectrum) - 10*(fr.ap_fit);
        spec_in = log10(nolog_eeg_psd(tmp_f,i,j));            
        fooof_diff = 10*(spec_in') - 10*(fr.ap_fit);
        tmp(:,i,j) = fooof_diff';
        fooof_freqs = fr.freqs;
    end
end
% figure;
% hold on;
% plot(log10(mean(squeeze(nolog_eeg_psd(:,:,j)),2)))
% plot(log10(mean(squeeze(nolog_eeg_psd(tmp_f,:,j)),2)));
% % plot(mean(squeeze(eeg_psd(:,:,j)),2))
% plot(fr.power_spectrum);
% plot(fr.ap_fit);
% plot(fr.fooofed_spectrum);
% hold off;
%-
eeg_psd = tmp;
freqs = fooof_freqs;
%## GET CLUSTER DATA
fprintf('Getting cluster information...\n');
tmp_subjs = {tmp_cl_study.datasetinfo.subject};
comp_arr = zeros(3,length(CLUSTER_PICS)); %(1, comps; 2, old_comps; 3, cluster
for i = 1:length(CLUSTER_PICS)
    cl_i = CLUSTER_PICS(i);
    ind = find(strcmp(EEG.subject,tmp_subjs));
    ind = tmp_study.cluster(cl_i).sets == ind;
    if any(ind)
        comp_arr(1,i) = tmp_study.cluster(cl_i).comps(ind);
        comp_arr(3,i) = cl_i;
        comp_arr(2,i) = EEG.etc.urreject.ic_keep(comp_arr(1,i));
    end
end
comp_arr = comp_arr(:,all(comp_arr,1));
%##
fprintf('Performing gait kinematic calculations...\n');
tband = (freqs > theta_band_lims(1) & freqs < theta_band_lims(2));
aband = (freqs > alpha_band_lims(1) & freqs < alpha_band_lims(2));
bband = (freqs > beta_band_lims(1) & freqs < beta_band_lims(2));
body_posx = find(contains({EEG.chanlocs.labels},'body_xpos','IgnoreCase',true));
body_posy = find(contains({EEG.chanlocs.labels},'body_ypos','IgnoreCase',true));
body_posz = find(contains({EEG.chanlocs.labels},'body_zpos','IgnoreCase',true));
EEG = eeg_checkset(EEG,'loaddata');
%##
steps_struct = repmat(def_step_structs,[1,length(EEG.epoch)*size(comp_arr,2)]);
cnt = 1;
%- loop through each condition
% errors: {'H1018'   }
% {'H1019'   }
% {'H1029'   }
% {'H1036'   }
% {'H1038'   }
% {'H1039'   }
% {'H2018_FU'}
% {'H2021'   }
% {'H2038'   }
% {'H2052'   }
% {'H3072'   }
% {'NH3070'  }
% {'NH3086'  }
% {'NH3090'  }
% {'NH3112'  }
% {'NH3113'  }
% {'NH3128'  }
for cond_i = 1:length(conds)
    %##
    inds = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{EEG.epoch.eventcond}));
    %- imu data
    datx = squeeze(EEG.data(body_posx,:,inds)); % AP (anteroposterior)
    daty = squeeze(EEG.data(body_posy,:,inds)); % ML (mediolateral)
    datz = squeeze(EEG.data(body_posz,:,inds)); % UD (Up-Down)
    %- sanity check
    %{
    figure;
    title(conds{cond_i})
    hold on;
    for s_i = 1:10
        plot3(datx(:,s_i),daty(:,s_i),datz(:,s_i));
    end
    hold off;
    %}
    %## ADVANCED EPOCH'ING?
    for s_i = 1:length(inds)-1
        tmp = EEG.epoch(inds(s_i));
        rhs_1 = nan();
        rhs_2 = nan();
        lhs_1 = nan();
        lto_1 = nan();
        rto_1 = nan();
        %## 
        t_rhi = find(strcmp(tmp.eventtype,'RHS'));
        lhi = find(strcmp(tmp.eventtype,'LHS'));
        rti = find(strcmp(tmp.eventtype,'RTO'));
        lti = find(strcmp(tmp.eventtype,'LTO'));
        %##
        % rht = sort(cell2mat(tmp.eventlatency(strcmp(tmp.eventtype,'RHS'))));
        % lht = sort(cell2mat(tmp.eventlatency(strcmp(tmp.eventtype,'LHS'))));
        % rtt = sort(cell2mat(tmp.eventlatency(strcmp(tmp.eventtype,'RTO'))));
        % ltt = sort(cell2mat(tmp.eventlatency(strcmp(tmp.eventtype,'LTO'))));
        %## DETERMINE LS FOOT EVENTS
        good_s = true;
        ii = 0;
        while good_s
            if length(t_rhi) > 2
                rhi = [t_rhi(1+ii),t_rhi(2+ii)];
            else
                rhi = t_rhi;
            end
            rht = [tmp.eventlatency{rhi}]*TIME_FACTOR;
            lht = [tmp.eventlatency{lhi}]*TIME_FACTOR;
            rtt = [tmp.eventlatency{rti}]*TIME_FACTOR;
            ltt = [tmp.eventlatency{lti}]*TIME_FACTOR;
            %-
            rhs_1 = rht(1); %tmp.eventlatency{rhi(1)}*TIME_FACTOR;
            rhs_2 = rht(2); %tmp.eventlatency{rhi(2)}*TIME_FACTOR;
            %- index critera
            % ind = lhi > rhi(1) & lhi < rhi(2);
            % lhs_1 = tmp.eventlatency{lhi(ind)}*TIME_FACTOR;
            % ind = lti > rhi(1) & lti < rhi(2);
            % lto_1 = tmp.eventlatency{lti(ind)}*TIME_FACTOR;
            % ind = rti > rhi(1) & rti < rhi(2);
            % rto_1 = tmp.eventlatency{rti(ind)}*TIME_FACTOR;
            %- latency criteria
            ind = lht > rht(1) & lht < rht(2);
            lhs_1 = lht(ind); %tmp.eventlatency{lhi(ind)}*TIME_FACTOR;
            ind = ltt > rht(1) & ltt < rht(2);
            lto_1 = ltt(ind); %tmp.eventlatency{lti(ind)}*TIME_FACTOR;
            ind = rtt > rht(1) & rtt < rht(2);
            rto_1 = rtt(ind); %tmp.eventlatency{rtt(ind)}*TIME_FACTOR;
            if all(~cellfun(@isempty,{rhs_1,rhs_2,lhs_1,lto_1,rto_1})) && length([rhs_1,rhs_2,lhs_1,lto_1,rto_1]) == 5
                %-
                log_s = true;
                good_s = false;
                fprintf('x');
            else
                % fprintf('Error. Trying next consecutive strikes...\n');
                fprintf('.');
                ii = ii + 1;
                if ii+2 > length(t_rhi)
                    fprintf('o');
                    good_s = false;
                    log_s = false;
                end
            end
        end
        %## CALCULATE IMU MEASURES
        %- rhs_1, lto_1, lhs_1, rto_1, rhs_2
        EVENT_ERR = 2.001; % potentially a millisecond discrepency for some reason
        if log_s
            %- rhs_1 ind
            diff = EEG.times - rhs_1/TIME_FACTOR;
            rhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
            %- rhs_1 ind
            diff = EEG.times - lto_1/TIME_FACTOR;
            ltoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
            %- rhs_1 ind
            diff = EEG.times - lhs_1/TIME_FACTOR;
            lhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
            %- rhs_1 ind
            diff = EEG.times - rto_1/TIME_FACTOR;
            rtoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
            %- rhs_2 ind
            diff = EEG.times - rhs_2/TIME_FACTOR;
            rhsi_2 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
            %- calculate pelvis movement
            %* step width (ml exc)
            vec = [daty(rhsi_1,s_i),daty(lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ml_exc = sqrt((m2-m1)^2);
            %* ap exc
            vec = [datx(rhsi_1,s_i),datx(lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ap_exc = sqrt((m2-m1)^2);
            %* ud exc
            vec = [datz(rhsi_1,s_i),datz(lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ud_exc = sqrt((m2-m1)^2);
            %* min max exc (ml)
            vec = [daty(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ml_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ap)
            vec = [datx(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ap_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ud)
            vec = [datz(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ud_mm_exc = sqrt((m2-m1)^2);
        end
        fprintf(' | ');
        %- struct
        if log_s %&& all(~isempty([rhs_1,rhs_2,lhs_1,lto_1,rto_1]))
            for i = 1:length(CLUSTER_PICS)
                cl_i = CLUSTER_PICS(i);
                comp_n = comp_arr(1,comp_arr(3,:) == cl_i);
                if ~isempty(comp_n)
                    steps_struct(cnt).subj_char = EEG.subject;
                    steps_struct(cnt).cond_char = conds{cond_i}; %speed_n_chars{cond_i}; %conds{cond_i};
                    steps_struct(cnt).group_char = EEG.group;
                    steps_struct(cnt).speed_char = speed_n_chars{cond_i}; %conds{cond_i};
                    steps_struct(cnt).speed_n = double(string(speed_n_chars{cond_i})); %conds{cond_i};
                    steps_struct(cnt).model_char = 'speed';
                    steps_struct(cnt).model_n = des_i;
                    steps_struct(cnt).cluster_n = cl_i;
                    steps_struct(cnt).comp_n = comp_n; % new comp number
                    steps_struct(cnt).comp_n_old = comp_arr(2,comp_arr(3,:) == cl_i); % old comp number
                    steps_struct(cnt).stride_n = inds(s_i);
                    steps_struct(cnt).rhs_1 = rhs_1;
                    steps_struct(cnt).rhs_2 = rhs_2;
                    steps_struct(cnt).lhs_1 = lhs_1;
                    steps_struct(cnt).lto_1 = lto_1;
                    steps_struct(cnt).rto_1 = rto_1;
                    steps_struct(cnt).gait_cycle_dur = rhs_2 - rhs_1;
                    steps_struct(cnt).stance_dur = rto_1 - rhs_1;
                    steps_struct(cnt).swing_dur = rhs_2 - rto_1;
                    steps_struct(cnt).double_sup_dur = rto_1 - lhs_1;
                    steps_struct(cnt).single_sup_dur = lhs_1 - lto_1;
                    steps_struct(cnt).step_dur = rhs_2 - lhs_1;
                    steps_struct(cnt).step_width = ml_exc;
                    steps_struct(cnt).ap_exc = ap_exc;
                    steps_struct(cnt).ud_exc = ud_exc;
                    steps_struct(cnt).step_width_mm = ml_mm_exc;
                    steps_struct(cnt).ap_exc_mm = ap_mm_exc;
                    steps_struct(cnt).ud_exc_mm = ud_mm_exc;
                    % step_structs(cnt).eeg_psd = squeeze(eeg_psd(:,cnt,comp_n));
                    % step_structs(cnt).freqs = freqs;
                    steps_struct(cnt).avg_theta = mean(eeg_psd(tband,inds(s_i),comp_n));
                    steps_struct(cnt).avg_alpha = mean(eeg_psd(aband,inds(s_i),comp_n));
                    steps_struct(cnt).avg_beta = mean(eeg_psd(bband,inds(s_i),comp_n));
                    steps_struct(cnt).avg_theta_post = mean(eeg_psd(tband,inds(s_i+1),comp_n));
                    steps_struct(cnt).avg_alpha_post = mean(eeg_psd(aband,inds(s_i+1),comp_n));
                    steps_struct(cnt).avg_beta_post = mean(eeg_psd(bband,inds(s_i+1),comp_n));
                    % step_structs(cnt).eeg_comp_arr = comp_arr;
                    % step_structs(cnt).eeg_comp_arr_deats = ['row 1 is the component in the reduced eeg set;' ...
                    %     ' row 2 is the component in the old eeg set; row 3 is the cluster number'];
                    cnt = cnt + 1;
                end
            end
        end
    end
    %- grab condition
    inds_cond = strcmp({steps_struct.cond_char},conds{cond_i});
    tmp_ss = steps_struct(inds_cond);
    %- per condition measures (step dur)
    musd = mean([tmp_ss.step_dur]);
    stdsd = std([tmp_ss.step_dur]);
    tmp = [tmp_ss.step_dur];
    %-
    vals = sqrt((tmp-musd).^2);
    vals = num2cell(vals);
    [steps_struct(inds_cond).var_step_dur_1] = deal(vals{:});
    %-
    vals = stdsd./sqrt((tmp-musd).^2);
    vals = num2cell(vals);
    [steps_struct(inds_cond).var_step_dur_2] = deal(vals{:});
    %-
    vals =sqrt((tmp-musd).^2)./stdsd;
    vals = num2cell(vals);
    [steps_struct(inds_cond).var_step_dur_3] = deal(vals{:});
    %- per condition measures (step width)
    musd = mean([tmp_ss.step_width_mm]);
    stdsd = std([tmp_ss.step_width_mm]);
    tmp = [tmp_ss.step_width_mm];
    %-
    vals = sqrt((tmp-musd).^2);
    vals = num2cell(vals);
    [steps_struct(inds_cond).var_step_width_1] = deal(vals{:});
    %-
    vals = stdsd./sqrt((tmp-musd).^2);
    vals = num2cell(vals);
    [steps_struct(inds_cond).var_step_width_2] = deal(vals{:});
    fprintf('\n');
end
%- assign to struct
% vals = num2cell(vals)';
% [steps_struct(inds_cond)] = deal(vals{:});
%- per design measures
%## SAVE    
steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.comp_n}));
steps_struct = struct2table(steps_struct);

end

