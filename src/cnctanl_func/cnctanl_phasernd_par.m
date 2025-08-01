function [C] = cnctanl_phasernd_par(EEG,conn_char,varargin)
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

%## PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'conn_char',@ischar)
%## PARAMETER
addOptional(p,'efm_cfg',def_efm_cfg)
addOptional(p,'cm_cfg',def_cm_cfg)
%##
parse(p,EEG,conn_char,varargin{:});
%## SET DEFAULTS
efm_cfg = p.Results.efm_cfg;
cm_cfg = p.Results.cm_cfg;
%% ===================================================================== %%
%## CREATE SURROGATE PHASE RAND
% Multiply each fourier amplitude by e^{iw)
% where w is a random phase chosen in [0 2pi]
% (c.f. Theiler, et al 1997)
% To generate hermitian phase distributions, we extract
% the phase of a random matrix. This ensures the surrogate
% spectrum is conjugate symmetric
EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
[npnts,nchs,ntr] = size(EEG.CAT.srcdata);
for tr=1:ntr
    EEG.CAT.srcdata(:,:,tr) = ...
        ifft(abs(fft(EEG.CAT.srcdata(:,:,tr))) ...
        .* exp(1i*angle(fft(rand(npnts,nchs)))), ...
        'symmetric');
end
EEG.CAT.srcdata = ipermute(EEG.CAT.srcdata,[2 1 3]);

%## Fit model and estimate connectivity for 1 permutation
%-- fit model
MODEL = est_fitMVAR('EEG',EEG,efm_cfg,'verb',0);

%-- calculate connectivity
cm_cfg.connmethods = {conn_char};
C = est_mvarConnectivity('EEG',EEG,'MODEL',MODEL,cm_cfg,'verb',0);

%## end
fprintf('done. cnctanl_phasernd_par.m: %0.2g min\n',toc(tt)/(60))
end