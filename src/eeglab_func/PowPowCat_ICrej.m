function bad_ics_out = PowPowCat_ICrej(EEG,varargin)
% ADAPTED BY CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% PowPowCat for IC rejection
% Reject ICs that have moderate cross-freq coupling in low (<8Hz) and high
% (>30Hz) frequency windows
% Takes the median correlation coefficient value in both high and low
% frequency windows, not including the identity coeffient (which always
% ==1)
%
% % Usage:
% >> [badPCC_IC] = PowPowCat_ICrej(EEG, varargin);
%
% Required inputs:
%   EEG            - EEG dataset structures with PowPowCat precomputed
%
% Optional inputs:
%   plotstuff      - plot extended component properties of "bad" ICs and
%                    spectral covariance (default = 0 (off))
%   outputFigureFolder - main folder path to store figures, subfolder will
%                       be created. (default = EEG.filepath)
% Output:
%   badPPC_IC      - Independent components identified as having moderate
%                    correlation in low and high freq. windows (e.g. eye
%                    blinks and muscle artifacts, respectively)
% Code Designers: Noelle Jacobsen, Jacob Salminen
% Code Date: 10/17/2021, MATLAB 2020b
% Authored (10/17/2021): Noelle Jacobsen, University of Florida, 
% edited (07/15/2023): Jacob Salminen, University of Florida; cleaned up 

%## TIME
tt = tic();
%## DEFINE DEFAULTS
%- Cutoff for Correlations ot determine bad components
% NOTE: correlation coefficient threshold, 0.5-0.7= moderate correlation
% REJ_STRUCT.ic_cutoff = 0.3; 
% REJ_STRUCT.low_freq_cutoff = 8;
% REJ_STRUCT.high_freq_cutoff = 30;
% NOTE (10/25/2023), Originally at 0.3; setting to 0.25 to be more
% aggressive on cleaning.
% NOTE (10/26/2023), Originally at 0.25; setting to 0.20 to be more
% aggressive on cleaning.
% NOTE (11/2/2023), Turning back to 0.3 for testing on new ICC
% parameterization.
% (11/7/2023), Turning to 0.1 & only performing on high frequency.
%-
%## Define Parser
p = inputParser;
DEF_REJ_STRUCT = struct('do_plot',false, ...
    'save_dir',{''}, ...
    'ic_cutoff',0.35, ...
    'low_freq_cutoff',8, ...
    'high_freq_cutoff',30);
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'REJ_STRUCT',DEF_REJ_STRUCT,@(x) validate_struct(x,DEF_REJ_STRUCT));
parse(p,EEG,varargin{:});
%## SET DEFAULTS
REJ_STRUCT = p.Results.REJ_STRUCT;
%-
REJ_STRUCT = set_defaults_struct(REJ_STRUCT,DEF_REJ_STRUCT);
%% ===================================================================== %%
%- extract powpowcat values
cov_matrix = EEG.etc.PowPowCAT.covMatrix;
freqs =  EEG.etc.PowPowCAT.freqs;
%- set params
lowfreq_coupling=cell(size(cov_matrix,3),1);
highfreq_coupling=cell(size(cov_matrix,3),1);
bad_ics_out=zeros(size(cov_matrix,3),1);
fprintf('\tIC#\t\tLowFreq Cov.\t\tHighFreq Cov.\n');
fprintf('________________________________________________');
for IC = 1:size(cov_matrix,3)
    %- extract correlations above and below alpha/beta frequencies (<8Hz,>30Hz)
    lowfreq_coupling{IC}=cov_matrix(freqs<REJ_STRUCT.low_freq_cutoff,freqs<REJ_STRUCT.low_freq_cutoff,IC); %below alpha
    highfreq_coupling{IC}=cov_matrix(freqs>REJ_STRUCT.high_freq_cutoff,freqs>REJ_STRUCT.high_freq_cutoff,IC); %above beta

    %## Low Freq Coupling
    CC = lowfreq_coupling{IC};
    %-- remove diagnoal identity (coeff always equals one with itself)
    CC(logical(eye(size(CC)))) =[]; 
    CC=reshape(CC,size(lowfreq_coupling{IC},1)-1,size(lowfreq_coupling{IC},1));
    %-- average (median)
    lowfreq_coupling_avg = median(median(CC));

    %## High Freq Coupling
    CC = highfreq_coupling{IC};
    %-- remove diagnoal identity (coeff always equals one with itself)
    CC(logical(eye(size(CC)))) =[]; 
    CC=reshape(CC,size(highfreq_coupling{IC},1)-1,size(highfreq_coupling{IC},1));
    %-- average (median)
    highfreq_coupling_avg = median(median(CC));
    
    %## TABLE
    fprintf('\n\t%i\t\t%.2f\t\t\t\t%.2f',IC,lowfreq_coupling_avg,highfreq_coupling_avg)
    %identify ICs with corr. coeff. above thresh in low and high freq windows
%     if lowfreq_coupling_avg>REJ_STRUCT.ic_cutoff || highfreq_coupling_avg>REJ_STRUCT.ic_cutoff
    if highfreq_coupling_avg > REJ_STRUCT.ic_cutoff
        bad_ics_out(IC) = IC;
        fprintf('\t**BAD**');
    end
end
bad_ics_out=bad_ics_out(bad_ics_out~=0);
%## TIME
fprintf('\ndone. PowPowCat_ICrej.m: %0.2g \n',toc(tt))
end

