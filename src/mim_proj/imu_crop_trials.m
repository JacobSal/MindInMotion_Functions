function [EEG_IMU] = imu_crop_trials(EEG_IMU, crop_xlsx_fpath, varargin)
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
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
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
%- logo
cat_logo();
%## TRIAL_CROPPING
cushion_sec = 5;
%% DEFINE DEFAULTS

%## PARSE
i = inputParser;
%## REQUIRED
addRequired(i,'EEG_IMU',@isstruct);
addRequired(i,'crop_xlsx_fpath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(i, EEG_IMU, crop_xlsx_fpath, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%% (CROP TRIAL DATA) =================================================== %%
subj_char = EEG_IMU.subject;
eeg_fname = EEG_IMU.filename;
tn = strsplit(eeg_fname,'.');
[ do_crop, crop_lats, do_crop_ls, crop_lats_ls ] = mim_check_trials(subj_char, ...
    tn{1},...
    crop_xlsx_fpath);

%-- remove gait events not in the ExactCrop_loadsol range
if do_crop_ls    
    % first convert the time --> latency
%             temp = [EEG_IMU.event.latency];
    fprintf('Cropping trial based on LoadSol data for %s...\n',tn{1});
    cushion_latency = EEG_IMU.srate*cushion_sec;
    crop_latency_range = EEG_IMU.srate*crop_lats_ls + cushion_latency + 1;
    total_gaitevent_idx = find(strcmp('LHS',{EEG_IMU.event.type}) | ...
        strcmp('LTO',{EEG_IMU.event.type})| ...
        strcmp('RHS',{EEG_IMU.event.type}) | ...
        strcmp('RTO',{EEG_IMU.event.type}));
    gait_event_keeps = cell(size(crop_latency_range,1),1);
    for i = 1:size(crop_latency_range,1)
        gait_event_keeps{i} = find([EEG_IMU.event.latency] > crop_latency_range(i,1) & ...
            [EEG_IMU.event.latency] < crop_latency_range(i,2));
    end
    gait_event_keeps = cat(2,gait_event_keeps{:});
    [gait_event_rmv] = setdiff(total_gaitevent_idx,gait_event_keeps);
    EEG_IMU.event(gait_event_rmv) = [];
end

%-- only deal with shortened trials for now (cannot ged rid of other events that happen during trial)
if do_crop 
    fprintf('Cropping trial based on researcher criteria %s...\n',tn{1});
    %--
    start_ind = find(strcmpi('TrialStart',{EEG_IMU.event.type}));
    % start_time = EEG_IMU.times(round(EEG_IMU.event(start_ind).latency))/1000;
    stop_ind = find(strcmpi('TrialEnd',{EEG_IMU.event.type}));
    crop_lats = EEG_IMU.srate*crop_lats(end)+EEG_IMU.event(start_ind).latency;        
    stop_time = EEG_IMU.times(round(crop_lats))/1000;
    %-- add trial end
    EEG_IMU.event(end+1) = EEG_IMU.event(stop_ind);
    EEG_IMU.event(end).latency = crop_lats;
    EEG_IMU.event(end).datetime = [];
    %--
    crop_update = crop_lats+cushion_sec;
    crop_update(1) = 0;
    crop_update(end,end) = stop_time+cushion_sec;
    EEG_IMU = pop_select( EEG_IMU, 'time',crop_update); 
    %(02/03/2025) ChangL? (JS), This is wrong in the original code from Roehl.
    % That code doesn't take into account the cushion time added for each trial
end