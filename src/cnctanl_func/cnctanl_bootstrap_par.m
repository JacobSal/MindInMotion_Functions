function [C] = cnctanl_bootstrap_par(EEG,resamp_idx,conn_char,varargin)
%CNCTANL_PHASERND_PAR this is an adaption from the SIFT function
%stat_surrogateGen
%gen
% Generate surrogate statistical distributions of SIFT connectivity and
% other estimators. This function can estimate bootstrap, jacknife, k-fold crossvalidation and
% phase-randomized (null) distributions.
%
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for 
%    updated version
% ===============================================
%
% === TODO ============================================
% - retest phaserand
% =====================================================
%
%
% Output                 Information                                                                                           
% -------------------------------------------------------------------------------------------------------------------------------
% PConn:                 Surrogate connectivity structure. Same format as EEG.CAT.Conn but with resamples stored in the last 
%                        dimension of each connectivity matrix (e.g. [M x M x Time x Freq x Resamples])   
%
% See Also: stat_analyticStats(), statcond(), stat_surrogateStats()
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 5.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu
%
% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% AMMENDED BY:::
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Jacob Salminen
% Code Date: 04/08/2025, MATLAB 2023b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%## TIME
tt = tic();
%## DEFINE DEFAULTS
%-- grab configs from MVAR fit
configs = EEG.CAT.configs;
def_efm_cfg = configs.est_fitMVAR;
def_cm_cfg = configs.est_mvarConnectivity;
def_bs_cfg = struct( ...
    'bs_type','continuos', ...
    'bs_cuts',1:10*EEG.srate:size(EEG.CAT.srcdata,2) ...
    );
%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'resamp_idx',@isnumeric);
addRequired(p,'conn_char',@ischar)
%## PARAMETER
addParameter(p,'bs_cfg',def_bs_cfg)
addParameter(p,'efm_cfg',def_efm_cfg)
addParameter(p,'cm_cfg',def_cm_cfg)
%##
parse(p,EEG,resamp_idx,conn_char,varargin{:});
%## SET DEFAULTS
bs_cfg = p.Results.bs_cfg;
efm_cfg = p.Results.efm_cfg;
cm_cfg = p.Results.cm_cfg;
%% ===================================================================== %%
%## CREATE SURROGATE BOOTTRAP 
% recalculates model several times and generates a distribution of fits
% all other bootstrapping/resampling modes

% select subset of trials
switch bs_cfg.bs_type
    case 'continuous'
        try
            tmp = cat(1,bs_cfg.bs_cuts{resamp_idx});
        catch e
            if strcmp(e.identifier,'MATLAB:catenate:dimensionMismatch')
                tmp = cat(2,bs_cfg.bs_cuts{resamp_idx});
            else
                error(getReport(e));
            end
        end
        EEG.CAT.srcdata = EEG.CAT.srcdata(:,tmp);
        EEG.CAT.trials = 1;
    otherwise
        EEG.CAT.srcdata = EEG.CAT.srcdata(:,:,resamp_idx);
        EEG.CAT.trials  = length(resamp_idx);
end

EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
EEG.CAT.srcdata = ipermute(EEG.CAT.srcdata,[2 1 3]);

%## Fit model and estimate connectivity for 1 permutation
%-- fit model
MODEL = est_fitMVAR('EEG',EEG,efm_cfg,'verb',0);

%-- calculate connectivity
cm_cfg.connmethods = {conn_char};
C = est_mvarConnectivity('EEG',EEG,'MODEL',MODEL,cm_cfg,'verb',0);

%## end
fprintf('done. cnctanl_bootstrap_par.m: %0.2g min\n',toc(tt)/(60))
end