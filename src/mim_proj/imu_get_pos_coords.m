function [EEG_IMU] = imu_get_pos_coords(EEG_IMU, crop_xlsx_fpath, varargin)
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
%## PATHS
if ~ispc
    pp = strsplit(path,':');
else
    pp = strsplit(path,';');
end
%## MATLAB FILTFILT
try
    %-
    inds = find(cellfun(@(x) contains(x,['toolbox' filesep 'signal' filesep 'signal']),pp));
    tmp_mat = dir(pp{inds(1)});
    fnames = {tmp_mat.name};
    inds_mat = cellfun(@(x) contains(x,'filtfilt'),fnames);
    fprintf('MATLAB filtfilt path: %s\n',[tmp_mat(inds_mat).folder filesep tmp_mat(inds_mat).name]);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('MATLAB filtfilt not found...\n');
    end
end
%## FIELDTRIP FILTFILT
try
    %-
    inds = find(cellfun(@(x) contains(x,['fieldtrip' filesep 'external' filesep 'signal']),pp));
    tmp_ft = dir(pp{inds(1)});
    fnames = {tmp_mat.name};
    inds_ft = cellfun(@(x) contains(x,'filtfilt'),fnames);
    fprintf('Fieldtrip filtfilt path: %s\n',[tmp_ft(inds_ft).folder filesep tmp_ft(inds_ft).name]);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('Fieldtrip filtfilt.m not found...\n');
    end
end
curr_ff = which('filtfilt');
fprintf('Current filtfilt path: %s\n',curr_ff);
%-
fprintf('Overriding using MATLAB filtfilt.m...\n')
try
    tmp = tmp_mat(inds).folder;
    path(tmp,path);
catch ME
    switch ME.identifier
        case 'MATLAB:needMoreRhsOutputs'
            fprintf('Path already set correctly for filtfilt.m .\n\n')
        otherwise
            fprintf('Error. %s\n%s\n\n',e.identifier,getReport(e));
    end
    
end
%% DEFINE DEFAULTS
if ~which('quaternRotate')
    error('You need to add the Gait-Tracking-With-x-IMU-master free library to get quaternion math');
end
%-- original parameters mim
% hpfilt_freq   = 0.2; 
hpfilt_freq =0.2;
%(01/24/2025) ryan (JS), Hi Colton. Ryan here. I decided to with 0.2 Hz instead of 0.25. Let me know how it looks
%(01/24/2025) JS, this may benefit from being upped to 0.3, but could also
%just be biasing the data too...
%(02/03/2025) JS, this may actually be over filtering. Could try 0.1 for
%now on
lpfilt_freq   = 20; %Hz
%-- updated?
% hpfilt_freq   = 0.15; %Hi Colton. Ryan here. I decided to with 0.2 Hz instead of 0.25. Let me know how it looks
% lpfilt_freq   = 10; %Hz
%(01/18/2025) JS, using parameterizations from: https://jneuroengrehab.biomedcentral.com/articles/10.1186/s12984-021-00828-0
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'EEG_IMU',@isstruct);
addRequired(p,'crop_xlsx_fpath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p, EEG_IMU, crop_xlsx_fpath, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%% CROP TRIAL DATA
EEG_STORE = EEG_IMU;

%## MANUAL CROP TRIALS??
% [EEG_IMU] = imu_crop_trials(EEG_IMU, crop_xlsx_fpath);
%(02/03/2025) JS, this doesn't seem to do much to improve the body position
%estimates (not noticeably by the eye).

%## TRIAL START & END DETERMINATION
start_i = find(strcmpi({EEG_IMU.event.type},'TrialStart'));
stop_i = find(strcmpi({EEG_IMU.event.type},'TrialEnd'));
start_lat = ceil(EEG_IMU.event(start_i).latency);
stop_lat = floor(EEG_IMU.event(stop_i).latency);
%-
% crop_times = [EEG_IMU.times(1), EEG_IMU.times(end)];
%-
crop_times = [EEG_IMU.times(start_lat)/1000, EEG_IMU.times(stop_lat)/1000];
EEG_IMU = pop_select(EEG_IMU,'time',crop_times);

%(02/03/2025) JS, pretty sure this is performed later on when
%the data is cropped using the laready cleaned ICA .set file.
%% Grab accelerometer channel data & APDM auto-generated quanternions
% regexp({EEG_IMU.chanlocs.labels},'acc');
%-
back_imu = contains({EEG_IMU.chanlocs.labels},'back','IgnoreCase',true);
acc_chan = contains({EEG_IMU.chanlocs.labels},'acc','IgnoreCase',true);
acc_chan = acc_chan & back_imu; %cellfun(@(x) ~isempty(x),acc_chan);
%- (optional) gyro, used to calculate orientation data
% gyro_chan = regexp({EEG_IMU.chanlocs.labels},'gyro'); 
% gyro_chan = cellfun(@(x) ~isempty(x),gyro_chan);
%- (optional) magnetometer, used to calculate orientation data
% mag_chan = regexp({EEG_IMU.chanlocs.labels},'mag'); 
% mag_chan = cellfun(@(x) ~isempty(x),mag_chan);
%-
orient_chan = contains({EEG_IMU.chanlocs.labels},'ori','IgnoreCase',true); 
orient_chan = orient_chan & back_imu; %cellfun(@(x) ~isempty(x),acc_chan);
%- converting to column data FYI
% sensor_frame_gyro = EEG_IMU.data(gyro_chan,:)'; % time_pntsx3 
% sensor_frame_mag = EEG_IMU.data(mag_chan,:)'; % time_pntsx3
sensor_frame_acc = EEG_IMU.data(acc_chan,:)'; % time_pntsx3 
%- Dataset containing the orientation quaternio
sensor_frame_ori = EEG_IMU.data(orient_chan,:)'; % time_pntsx4
%- remove nan's
sensor_frame_acc(isnan(sensor_frame_acc)) = 0;
sensor_frame_ori(isnan(sensor_frame_ori)) = 0;
fprintf('percent accelerometer data nan: %0.2g\n',sum(logical(isnan(sensor_frame_acc)),'all')/(size(sensor_frame_acc,1)*size(sensor_frame_acc,2)));
fprintf('percent orientation data nan: %0.2g\n',sum(logical(isnan(sensor_frame_ori)),'all')/(size(sensor_frame_ori,1)*size(sensor_frame_ori,2)));
%% GET SENSOR POSITION DATA
sensor_frame_acc = EEG_IMU.data(acc_chan,:)'; % time_pntsx3 
sensor_frame_acc_nograv = sensor_frame_acc;
%- de-mean sensor
sensor_frame_acc_nograv = sensor_frame_acc_nograv - mean(sensor_frame_acc_nograv);
% sensor_frame_acc_nograv(:,1) = sensor_frame_acc_nograv(:,1) - mean(sensor_frame_acc_nograv(:,1));
%- (OPTIONAL) filter accelerometer data
sensor_frame_acc_nograv = FilterIMUfunction(sensor_frame_acc_nograv,0.3,20, EEG_IMU.srate);
% sensor_frame_acc_nograv(:,3) = FilterIMUfunction(sensor_frame_acc_nograv(:,3),0.15,10, EEG_IMU.srate);
%- (OPTIONAL) 
% dig_filt = designfilt('highpassiir', ...
%     'FilterOrder',4,...
%     'HalfPowerFrequency',0.3, ...
%     'DesignMethod','butter', ...
%     'SampleRate',EEG_IMU.srate); % high pass iir
% sensor_frame_acc_nograv = filtfilt(dig_filt,sensor_frame_acc_nograv);
%(01/14/2025) JS, this may cause issues if you filter acc, vel, and pos.
%- Integrate acceleration to obtain velocity
sensor_frame_vel = cumtrapz((sensor_frame_acc_nograv))*(1/EEG_IMU.srate); %m/s
%- filter velocity data
sensor_frame_vel = FilterIMUfunction(sensor_frame_vel,hpfilt_freq,lpfilt_freq, EEG_IMU.srate);
%(3/24/2021) colton (JS), this seems to improve position generation
%validity. Not necessarily accuracy.
%(03/24/2021) RD (JS), may not be needed to detrend given we are now BP
%filtering the velocity.
%- integrate velocity to obtain world position
sensor_frame_pos = (cumtrapz((sensor_frame_vel))*(1/EEG_IMU.srate) ); %m
%- filter position data
sensor_frame_pos = FilterIMUfunction(sensor_frame_pos,hpfilt_freq,lpfilt_freq, EEG_IMU.srate);

%## SANITY CHECK
%- state plots
% index 1 = up-down
% index 2 = right-left
% index 3 = ant-post
fig = figure;
hold on;
subplot(1,3,1); 
plot(sensor_frame_pos(:,2),sensor_frame_pos(:,3),':'); 
xlabel('Right'); 
ylabel('Foward'); 
title('Transverse'); %top down 
subplot(1,3,2); 
plot(sensor_frame_pos(:,3),sensor_frame_pos(:,1),':'); 
xlabel('Forward'); 
ylabel('Up'); 
title('Sagittal'); %side view
subplot(1,3,3); 
plot(sensor_frame_pos(:,2),-sensor_frame_pos(:,1),':'); 
xlabel('Right'); 
ylabel('Up'); 
title('Coronal'); %from behind or in front
sgtitle('sensor positions')
hold off;
nn = strsplit(EEG_IMU.filename,'.');
% saveas(fig,[EEG_IMU.filepath filesep sprintf('sensor_position_%s.fig',nn{1})])
% close(fig)
%{
%- 2d timeseries
imu_times = EEG_IMU.times/1000;
fig = plot_imu_signals(sensor_frame_acc_nograv,sensor_frame_vel,sensor_frame_pos,imu_times)
%- 3d
figure;
datx = sensor_frame_pos(:,1);
daty = sensor_frame_pos(:,2); 
datz = sensor_frame_pos(:,3); 
hold on;
p = plot3(datx(:,:),daty(:,:),datz(:,:),'LineStyle',':');
p.Color = [p.Color, 0.75];
view([45,45]);
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
%}
%% CONVERT SENSOR ACC. DATA TO WORLD POSITION DATA
% hpfilt_freq =0.1;
%## APPLY OPAL QUANTERNIONS TO SENSOR ACC DATA
sensor_frame_acc = EEG_IMU.data(acc_chan,:)'; % time_pntsx3 
%-- (OPTIONAL) de-mean sensor
% sensor_frame_acc = sensor_frame_acc - mean(sensor_frame_acc);
%- (OPTIONAL) filter accelerometer data (remove drifts)
% sensor_frame_acc = FilterIMUfunction(sensor_frame_acc,0.3,20, EEG_IMU.srate);
% sensor_frame_acc(:,3) = FilterIMUfunction(sensor_frame_acc(:,3),hpfilt_freq,lpfilt_freq, EEG_IMU.srate);
%- (OPTIONAL) remove acc. drifts *** this option could help, but also
%seemingly unnecessary, no literature to support
% acc_means = mean(sensor_frame_acc);
% sensor_frame_acc = sensor_frame_acc - acc_means;
% dig_filt = designfilt('highpassiir', ...
%     'FilterOrder',4,...
%     'HalfPowerFrequency',0.03, ...
%     'DesignMethod','butter', ...
%     'SampleRate',EEG_IMU.srate); % high pass iir
% sensor_frame_acc = filtfilt(dig_filt,sensor_frame_acc);
% sensor_frame_acc = sensor_frame_acc + acc_means;
%-- (OPTION) 
% world_frame_acc = zeros(size(sensor_frame_acc));
% window_sz = 1;
% sz = 1;
% while sz < size(sensor_frame_acc,1)-(window_sz+1)
%     quants = mean(sensor_frame_ori(sz:sz+window_sz,:));
%     quants = repmat(quants,[window_sz,1]);
%     world_frame_acc(sz:sz+window_sz-1,:) = quaternRotate(sensor_frame_acc(sz:sz+window_sz-1,:), quants);
%     sz = sz + window_sz;
%     % disp(sz);
% end
% sz_dif = size(sensor_frame_acc,1)-sz;
% quants = mean(sensor_frame_ori(sz:sz+sz_dif,:));
% quants = repmat(quants,[sz_dif,1]);
% world_frame_acc(sz:sz+sz_dif-1,:) = quaternRotate(sensor_frame_acc(sz:sz+sz_dif-1,:), quants);
%-- (OPTION) go through each time point
% world_frame_acc = zeros(size(sensor_frame_acc));
% for samp_i = 1:size(sensor_frame_acc,1)
%     %- grab the quaternion that tells us how our local IMU frame relates to the world frame (sensor fusion already done)
%     quants = sensor_frame_ori(samp_i,:); 
%     %- rotate the local accelerations to world frame
%     world_frame_acc(samp_i,:)= quaternRotate(sensor_frame_acc(samp_i,:), quants); 
% end
%-- (OPTION) use vector math to perform transformation
world_frame_acc = quaternRotate(sensor_frame_acc, sensor_frame_ori);
%(01/14/2025) JS, not sure why the previous method did sample by sample,
%unsure if it matters. This way would save comp time.

%## REMOVE NAN'S BEFORE FILTERING
%(12/01/2024) JS, bug where filtering will not handle NaN's well and will
%NaN whole trace if they are not removed.
world_frame_acc(isnan(world_frame_acc)) = deal(0);
fprintf('percent data points nan: %0.2g\n',sum(logical(isnan(world_frame_acc)),'all')/(size(world_frame_acc,1)*size(world_frame_acc,2)));

%## REMOVE GRAVITY EFFECTS
g_const = mean(world_frame_acc(:,3)); %units m/s^2 %should be positive or negative????
fprintf('Gravity = %s m/s^2\n',num2str(g_const));
%-- demean accelerometer data?
world_frame_acc_nograv = world_frame_acc - mean(world_frame_acc);
%-- Optionally filter the accelerometer
% world_frame_acc_nograv = FilterIMUfunction(world_frame_acc_nograv,hpfilt_freq,lpfilt_freq, EEG_IMU.srate);
%(01/14/2025) JS, this may cause issues if you filter acc, vel, and pos.
%-- Integrate acceleration to obtain velocity
world_frame_vel = cumtrapz((world_frame_acc_nograv))*(1/EEG_IMU.srate); %m/s
%-- filter velocity data
world_frame_vel = FilterIMUfunction(world_frame_vel,hpfilt_freq,lpfilt_freq, EEG_IMU.srate);
%(3/24/2021) colton (JS), this seems to improve position generation
%validity. Not necessarily accuracy.
%(03/24/2021) RD (JS), may not be needed to detrend given we are now BP
%filtering the velocity.
%-- integrate velocity to obtain world position
world_frame_pos = (cumtrapz((world_frame_vel))*(1/EEG_IMU.srate) ); %m
%-- filter position data
world_frame_pos = FilterIMUfunction(world_frame_pos,hpfilt_freq,lpfilt_freq, EEG_IMU.srate);

%## SANITY CHECK
%- state plots
fig = figure;
hold on;
subplot(1,3,1); 
plot(-world_frame_pos(:,2),world_frame_pos(:,1),':'); 
xlabel('Right'); ylabel('Foward'); 
title('Transverse'); %top down 
subplot(1,3,2); 
plot(world_frame_pos(:,1),world_frame_pos(:,3),':'); 
xlabel('Forward'); 
ylabel('Up'); 
title('Sagittal'); %side view
subplot(1,3,3); 
plot(-world_frame_pos(:,2),world_frame_pos(:,3),':'); 
xlabel('Right'); 
ylabel('Up'); 
title('Coronal'); %from behind or in front
sgtitle('world positions')
hold off;
nn = strsplit(EEG_IMU.filename,'.');
saveas(fig,[EEG_IMU.filepath filesep sprintf('world_position_%s.fig',nn{1})])
% close(fig)
%{
%- 2d timeseries
imu_times = EEG_IMU.times/1000;
fig = plot_imu_signals(world_frame_acc_nograv,world_frame_vel,world_frame_pos,imu_times)
%- 3d
figure;
datx = world_frame_pos(:,1);
daty = world_frame_pos(:,2); 
datz = world_frame_pos(:,3); 
hold on;
p = plot3(datx(:,:),daty(:,:),datz(:,:),'LineStyle',':');
p.Color = [p.Color, 0.75];
view([45,45]);
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
%- state plots
fig_i = figure;
hold on;
subplot(1,3,1); 
plot(-world_frame_pos(:,2),world_frame_pos(:,1),':'); 
xlabel('Right'); ylabel('Foward'); 
title('Transverse'); %top down 
subplot(1,3,2); 
plot(world_frame_pos(:,1),world_frame_pos(:,3),':'); 
xlabel('Forward'); 
ylabel('Up'); 
title('Sagittal'); %side view
subplot(1,3,3); 
plot(-world_frame_pos(:,2),world_frame_pos(:,3),':'); 
xlabel('Right'); 
ylabel('Up'); 
title('Coronal'); %from behind or in front
hold off;
%}
%% GET BODY FRAME POSITIONS
% hpfilt_freq =0.15;
%-- find avg body frame directions from sensor orientations (front/left/up)
body_to_world_trans = mean(sensor_frame_ori); % go from sensor to body(world?) frame
%(01/18/2025) JS, the original code does not have the most consistent
%documentaiton, so trying to fix naming conventions and comments as I go.
%This pariticular section of code is very confusing. Okay, typically body
%position of the pelvis should be almost 1-to-1 with the sensor, but
%treadmill data makes deriving forward directino difficult, so by taking
%you can generally define the forward direction of the body as the
%"strongest" forward movement vector of the orientation data. This is the
%thought, but could be invalid in some scenarios. 

%-- get world to body frame filtered data?
% world_to_body_trans = quaternConj(body_to_world_trans);
%-- get body to world frame transformation matrix
body_to_world_fwd           = quaternRotate([0 0 -1],body_to_world_trans);
body_to_world_left          = quaternRotate([0 -1 0],body_to_world_trans);
% body_to_world_grav          = quaternRotate([-1 0 0],body_to_world_trans);

%-- original rotmat calculation
body_to_world_rotmat        = [[(body_to_world_fwd(1:2)/sqrt(sum(body_to_world_fwd(1:2).^2))),0]; ...
    [(body_to_world_left(1:2)/sqrt(sum(body_to_world_left(1:2).^2))),0]; ...
    [0,0,1]];
%(01/25/2025) JS, the original line: 
%       (body_to_world_fwd(1:2)/sqrt(sum(body_to_world_fwd(1:2).^2)))
% does not use a element wise division as oppose to the second line:
%       (body_to_world_left(1:2)./sqrt(sum(body_to_world_left(1:2).^2)))
% unsure if this is intentional or not. trying out things... the
% denominator is singular. the element wise vs. matrix division doesn't
% matter

%-- try 2 rotmat calc.
% body_to_world_rotmat        = [[(body_to_world_fwd(1:2)/sqrt(sum(body_to_world_fwd(1:2).^2))),0]; ...
%     [(body_to_world_left(1:2)./sqrt(sum(body_to_world_left(1:2).^2))),0]; ...
%     body_to_world_grav];

%-- a more sensible rotmat calc?
% body_to_world_rotmat = [body_to_world_fwd; body_to_world_left; body_to_world_grav];
%(01/18/2025) JS, some weird code from old function: 
%2x2 %pca rotation matrix. score = centeredData*coeff; centeredData = SCORE*COEFF'
% zRotAngle = acos(body_to_world_rotmat(1,1)); %analytic solution
% tempQ = rotMat2quatern( euler2rotMat(0,0,-zRotAngle) ); 
% negative to flip reference frame?  %Needs Gait-Tracking-.... library and subfolders to be added to path

%- find world to body transformation from the inverse of the world to body
world_to_body_rotmat = pinv(body_to_world_rotmat);
%(03/16/2021) ryan (JS), ryan update in future
%(01/18/2025) JS, this did not get updated... It simply flips the signs of
%the diagonals

%- rotate all data to avg body frame of reference
pos_est_body = world_frame_pos*world_to_body_rotmat; %Front, Left, Up
vel_est_body = world_frame_vel*world_to_body_rotmat; %Front, Left, Up
acc_est_body = world_frame_acc_nograv*world_to_body_rotmat; %Front, Left, Up
%(01/18/2025) ryan (JS), rotate things from global frame to something more like A/P and mediolateral

%## SANITY CHECK
%- state plots
fig = figure;
hold on;
subplot(1,3,1); 
plot(-pos_est_body(:,2),pos_est_body(:,1),':'); 
xlabel('Right'); ylabel('Foward'); 
title('Transverse'); %top down 
subplot(1,3,2); 
plot(pos_est_body(:,1),pos_est_body(:,3),':'); 
xlabel('Forward'); 
ylabel('Up'); 
title('Sagittal'); %side view
subplot(1,3,3); 
plot(-pos_est_body(:,2),pos_est_body(:,3),':'); 
xlabel('Right'); 
ylabel('Up'); 
title('Coronal'); %from behind or in front
hold off;
sgtitle('body positions')
nn = strsplit(EEG_IMU.filename,'.');
% saveas(fig,[EEG_IMU.filepath filesep sprintf('body_position_%s.fig',nn{1})])
% close(fig)
%{
%- 2d timeseries
imu_times = EEG_IMU.times/1000;
fig = plot_imu_signals(acc_est_body,vel_est_body,pos_est_body,imu_times)
%- 3d
figure;
datx = pos_est_body(:,1);
daty = pos_est_body(:,2); 
datz = pos_est_body(:,3); 
hold on;
p = plot3(datx(:,:),daty(:,:),datz(:,:),'LineStyle',':');
p.Color = [p.Color, 0.75];
view([45,45]);
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
%}
%% make static image of entire timecourse of trial from 3 orthogonal angles
%{
fig_i = figure;
hold on;
subplot(1,3,1); 
plot(-pos_est_body(:,2),pos_est_body(:,1),':'); 
xlabel('Right'); ylabel('Foward'); 
title('Transverse'); %top down 
subplot(1,3,2); 
plot(pos_est_body(:,1),pos_est_body(:,3),':'); 
xlabel('Forward'); 
ylabel('Up'); 
title('Sagittal'); %side view
subplot(1,3,3); 
plot(-pos_est_body(:,2),pos_est_body(:,3),':'); 
xlabel('Right'); 
ylabel('Up'); 
title('Coronal'); %from behind or in front
hold off;
saveas(fig_i,[sprintf('state_space_flu_%s.jpg',tname)]);
%}
%% IMU 6DOF MOVIE (takes a lot of time to render)
%{
span = 10000; % display 10000 ms of data
s_i = randi(size(pos_est_body,1)-(span+1),1);
s_i = s_i:s_i+span;
imu_times = EEG.times(s_i);
%- prep
pos_in = pos_est_body(s_i,:);
rot_in = quatern2rotMat(quaternConj(sensor_frame_ori(s_i,:)));
disp_imu_movie(pos_in,rot_in,[save_dir filesep sprintf('6DOF_movie_%s',trial_name)]);
%}
%% MAKE NEW EEG SET WITH NEW BODY POSITION DATA
%## BUFFER NEW CHANS
%(12/13/2024) JS, altering to make sure that start and end of trial isn't
%messing up the overall alignment of the IMU (simply recropping in the data
%after alignment with the cropped data)
inds_crop = EEG_STORE.times >= crop_times(1)*1000 & EEG_STORE.times <= crop_times(2)*1000;
diff_vals = size(EEG_IMU.data,2) - sum(inds_crop);
intv = 1;
cnt = 0;
while diff_vals ~= 0 && cnt < 100
    add_inds = find(inds_crop == 0);
    diffdiff = find(diff(add_inds)>1);
    inds_crop(diffdiff(intv)) = true;
    diff_vals = size(EEG_IMU.data,2) - sum(inds_crop);
    if intv == 1 && length(diffdiff) > 1
        intv = 2;
    else
        intv = 1;
    end
    cnt = cnt + 1;
end

%## ASSIGNING NEW DATA
tmp_data = zeros(size(EEG_STORE.data,1)+16,size(EEG_STORE.data,2));
n_chans = size(EEG_STORE.data,1);
tmp_data(1:n_chans,inds_crop) = EEG_IMU.data;
tmp_data(n_chans+1:n_chans+3,inds_crop) = pos_est_body';
tmp_data(n_chans+4:n_chans+6,inds_crop) = vel_est_body';
tmp_data(n_chans+7:n_chans+9,inds_crop) = acc_est_body';
tmp_data(n_chans+10:n_chans+12,inds_crop) = world_frame_pos';
tmp_data(n_chans+13:n_chans+15,inds_crop) = sensor_frame_pos';
tmp_data(n_chans+16:n_chans+19,inds_crop) = sensor_frame_ori';

%## EXPORT EEG
EEG_IMU = EEG_STORE;
n_chans = size(EEG_STORE.data,1);
EEG_IMU.pnts = size(EEG_IMU.data,2);
EEG_IMU.etc.body_to_world_trans = body_to_world_trans;
EEG_IMU.etc.world_to_body_trans = world_to_body_rotmat;
EEG_IMU.data = tmp_data;
EEG_IMU.chanlocs(n_chans+1).labels  = 'body_xpos';
EEG_IMU.chanlocs(n_chans+2).labels  = 'body_ypos';
EEG_IMU.chanlocs(n_chans+3).labels  = 'body_zpos';
EEG_IMU.chanlocs(n_chans+4).labels  = 'body_xVel';
EEG_IMU.chanlocs(n_chans+5).labels  = 'body_yVel';
EEG_IMU.chanlocs(n_chans+6).labels  = 'body_zVel';
EEG_IMU.chanlocs(n_chans+7).labels  = 'body_xAcc';
EEG_IMU.chanlocs(n_chans+8).labels  = 'body_yAcc';
EEG_IMU.chanlocs(n_chans+9).labels  = 'body_zAcc';
EEG_IMU.chanlocs(n_chans+10).labels  = 'world_xpos';
EEG_IMU.chanlocs(n_chans+11).labels  = 'world_ypos';
EEG_IMU.chanlocs(n_chans+12).labels  = 'world_zpos';
EEG_IMU.chanlocs(n_chans+13).labels  = 'sensor_xpos';
EEG_IMU.chanlocs(n_chans+14).labels  = 'sensor_ypos';
EEG_IMU.chanlocs(n_chans+15).labels  = 'sensor_zpos';
EEG_IMU.chanlocs(n_chans+16).labels = 'q1';
EEG_IMU.chanlocs(n_chans+17).labels = 'qi';
EEG_IMU.chanlocs(n_chans+18).labels = 'qj';
EEG_IMU.chanlocs(n_chans+19).labels = 'qk';
for i = 1:length(EEG_IMU.chanlocs)
    if contains(EEG_IMU.chanlocs(i).labels,'Back')
        out = strsplit(EEG_IMU.chanlocs(i).labels,'_');
        out = strjoin(['orig',out(2:end)],'_');
        EEG_IMU.chanlocs(i).labels = out;
    end
    % NEW_EEG_IMU.chanlocs(i).type = 'IMU_Opal';
    EEG_IMU.chanlocs(i).type = 'BioM';
end
EEG_IMU.urchanlocs = EEG_IMU.chanlocs;
EEG_IMU.nbchan = length(EEG_IMU.chanlocs);
%##
toc
end

%% SUBFUNCTIONS
function [dat] = FilterIMUfunction(dat,hpfilt_freq,lpfilt_freq, srate)
    %This functino expects column vectors
    %It applies 4th order butterworth filtering (forward and reverse with filtfilt)
    %You can do BP filtering with this function. 
    %You can also do HP only if you set the LP input to inf.
    %This function does some special version of padding prior to filtering to prevent edge effects 
    
    %## filter position data
    %note: we are doing fancy copying of fake data to the left and right of
    %original dataset so we can filter the data without edge effects. Our
    %original IMU data is sharply cut into 180-sec long trials. We want to use
    %all of that data, but we also need to filter. Since filters use windows,
    %we need extra padded data.
    
    %temporarily add extra data to edges to help avoid edge effects
    %borrowing code found online for minimizing edge effects
    R=0.1; % 10% of signal
    Nr=5000;
    N=size(dat,1);
    % At most xxx points
    NR=min(round(N*R),Nr);
    %--
    x1 = zeros(length(2:NR+1),size(dat,2));
    x2 = zeros(length(2:NR+1),size(dat,2));
    %-- maintain continuity in level and slope
    for i=1:size(dat,2)
        x1(:,i) = 2*dat(1,i) - flipud(dat(2:NR+1,i));  
        x2(:,i) = 2*dat(end,i) - flipud(dat(end-NR:end-1,i));
    end
    %-- append extra data to each size of the original (x)
    dat=[x1;dat;x2]; 

    %## DESIGN FILTER
    if isinf(lpfilt_freq) % highpass only
        fprintf('Butterworth highpass filtering\n');
        fprintf('[%0.2f,Inf] Hz\n',hpfilt_freq);
        % dig_filt = designfilt('highpassiir', 'FilterOrder', 4,...
        %                         'HalfPowerFrequency', hpfilt_freq, ...
        %                         'SampleRate',srate); % hp only
        dig_filt = designfilt('highpassiir', 'FilterOrder', 4,...
                                'HalfPowerFrequency', hpfilt_freq, ...
                                'DesignMethod','butter', ...
                                'SampleRate',srate); % hp only
    else % if HP and LP cutoffs provided (non-inf)
        fprintf('Butterworth bandpass filtering...\n');
        fprintf('[%0.2f,%0.2f] Hz\n',hpfilt_freq, lpfilt_freq);
        % dig_filt = designfilt('bandpassiir', 'FilterOrder', 4,...
        %                         'HalfPowerFrequency1', hpfilt_freq, ...
        %                         'HalfPowerFrequency2', lpfilt_freq, ...
        %                         'SampleRate', srate); % bandpass
        dig_filt = designfilt('bandpassiir', 'FilterOrder', 4,...
                                'HalfPowerFrequency1', hpfilt_freq, ...
                                'HalfPowerFrequency2', lpfilt_freq, ...
                                'DesignMethod','butter', ...
                                'SampleRate', srate); % bandpass
    end
    
    %## APPLY FILTER    
    %- implementation with updated design filt function adapted for fieldtrip 
    % tmp = filtfilt(dig_filt.Coefficients(1,:),dig_filt.Coefficients(2,:),dat); 
    dat = filtfilt(dig_filt,dat);
    %-
    % tmpf = filter(dig_filt,dat); % my update implement?
    % tmp = filtfilt(dig_filt,dat); % ryan's/natalie's implementation
    %(02/09/2023) JS, changed to using filter as recommended by the designfilt function
    %(01/15/2025) JS, need to use filtfilt as it is a zero-phase filtering
    %method that prevents shifts in time after filtering.
    %- crop back to original data (temporarily added extra to the front and back)
    dat = dat(NR+1:end-NR,:);
end

function fig = plot_imu_signals(acc,vel,pos,imu_times)
    y_shift = 0.25;
    ax_y = 0.75;
    ax_x = 0.1;
    ax_h = 0.2;
    ax_w = 0.75;
    AXES_DEFAULT_PROPS = {'box','off', ...
        'xtick',[], ...
        'ytick',[], ...
        'ztick',[], ...
        'xcolor',[1,1,1], ...
        'ycolor',[1,1,1]};
    labs = {'x','y','z'};
    LINE_ALPHA = 0.75;
    pp = [];
    %-
    fig = figure('color','white');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    %-
    ax = axes();
    hold on;
    for i = 1:size(acc,2)
        p = plot(ax,imu_times,acc(:,i),'LineStyle',':','LineWidth',1,'DisplayName',sprintf('%s',labs{i}));
        p.Color = [p.Color, LINE_ALPHA];
        pp = [pp,p];
    end
    set(ax,'Position',[ax_x,ax_y,ax_w,ax_h]);  %[left bottom width height]
    ylabel(ax,'Acceleration (m/s^2)');
    xlabel(ax,'Time (s)');    
    %-
    ax = axes();
    hold on;
    for i = 1:size(vel,2)
        p = plot(ax,imu_times,vel(:,i),'LineStyle',':','LineWidth',1,'DisplayName',sprintf('%s',labs{i}));
        p.Color = [p.Color, LINE_ALPHA];
    end
    set(ax,'Position',[ax_x,ax_y-y_shift,ax_w,ax_h]);  %[left bottom width height]
    ylabel(ax,'Velocity (m/s)');
    xlabel(ax,'time (s)');
    %-
    ax = axes();
    hold on;
    for i = 1:size(pos,2)
        p = plot(ax,imu_times,pos(:,i),'LineStyle',':','LineWidth',1,'DisplayName',sprintf('%s',labs{i}));
        p.Color = [p.Color, LINE_ALPHA];
    end
    set(ax,'Position',[ax_x,ax_y-y_shift*2,ax_w,ax_h]);  %[left bottom width height]
    ylabel(ax,'Position (m)');
    xlabel(ax,'time (s)');
    legend(pp);
    hold off;

end