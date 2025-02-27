% Modified by Chang Liu - 20220714
% Crop trials based on 'TrialCrop' sheet
% Remove gait events based on 'Loadsol' sheet
function [ DoCrop, ExactCrop, DoCrop_loadsol, ExactCrop_loadsol ] = CropTrialCheckFunc_checkLoadsol(SubjStr,TrialName,inputdirectory)
    %IgnoreTrialCheckFunc - looks through a master Excel sheet to see if we
    %need to ignore the current trial
    %   Detailed explanation goes here
    if isempty(inputdirectory)
        disp('couldn''t find the excel doc so the function can''t work')
    end
    
    MasterTable = readtable(inputdirectory,'Sheet','trialCrop','Range','B2:V127'); %loads it up
    MasterTable.Properties.VariableNames(1:21) = {'sorting_num','SubjCode','SP_0p5_1',...
        'SP_0p5_2',	'SP_0p25_1','SP_0p25_2','SP_0p75_1','SP_0p75_2','SP_1p0_1',	'SP_1p0_2',	...
        'TM_flat_1','TM_flat_2','TM_high_1','TM_high_2','TM_low_1',	'TM_low_2',	'TM_med_1',	'TM_med_2',	'Rest','MotorImagery_1','MotorImagery_2'};

    SubjectMatch = strcmp(MasterTable{:,'SubjCode'},SubjStr); %look thru the subj to find match
    SubjectMatchInd = find(SubjectMatch); %grab corresponding index

    if isempty(SubjectMatchInd)
        DoCrop = false; %haven't ran the trial so there is not even a match
        ExactCrop = []; %2021-12-07 RJD added to avoid error in analyzeLoadsolManySubj script
    else
        Result = MasterTable{SubjectMatchInd,TrialName}; %result as a cell
        ResultString = Result{1}; %the result as a string
        ExactCrop = str2num(ResultString); %result as a vector. the exact seconds to crop it by
        %should get something in format of [boundary boundary]
        
        %2021-06-16 RJD replaced DoCrop  logic because it could not handle
        %excel having an empty cell (i.e. ResultString = '', ExactCrop = [])
        if isempty(ExactCrop)
            DoCrop = false;
        else
            DoCrop = true;
        end
            
%         if ResultString(1) == '[' && ResultString(length(ResultString)) == ']' %kinda bad logic but also kinda good 
%             DoCrop = true;
%         else
%             DoCrop = false;
%         end
    
    end
    
    %% 
    LoadsolTable = readtable(inputdirectory,'Sheet','loadsolCrop'); %loads it up
    LoadsolTable.Properties.VariableNames(1:21) = {'sorting_num','SubjCode','SP_0p5_1',...
        'SP_0p5_2',	'SP_0p25_1','SP_0p25_2','SP_0p75_1','SP_0p75_2','SP_1p0_1',	'SP_1p0_2',	...
        'TM_flat_1','TM_flat_2','TM_high_1','TM_high_2','TM_low_1',	'TM_low_2',	'TM_med_1',	'TM_med_2',	'Rest','MotorImagery_1','MotorImagery_2'};
    SubjectMatch = strcmp(LoadsolTable{:,'SubjCode'},SubjStr); %look thru the subj to find match
    SubjectMatchInd = find(SubjectMatch); %grab corresponding index

    if isempty(SubjectMatchInd)
        DoCrop_loadsol = false; %haven't ran the trial so there is not even a match
        ExactCrop_loadsol = []; %2021-12-07 RJD added to avoid error in analyzeLoadsolManySubj script
    else
        Result = LoadsolTable{SubjectMatchInd,TrialName}; %result as a cell
        ResultString = Result{1}; %the result as a string
        ExactCrop_loadsol = str2num(ResultString); %result as a vector. the exact seconds to crop it by
        %should get something in format of [boundary boundary]
        
        %2021-06-16 RJD replaced DoCrop  logic because it could not handle
        %excel having an empty cell (i.e. ResultString = '', ExactCrop = [])
        if isempty(ExactCrop_loadsol)
            DoCrop_loadsol = false;
        else
            DoCrop_loadsol = true;
        end
            
%         if ResultString(1) == '[' && ResultString(length(ResultString)) == ']' %kinda bad logic but also kinda good 
%             DoCrop = true;
%         else
%             DoCrop = false;
%         end
    
    end
end

