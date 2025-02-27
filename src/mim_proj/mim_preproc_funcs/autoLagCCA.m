function [EEG_clean] = autoLagCCA(EEG,lagAmount_samples,Rsq_thres)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

EEG_clean = EEG;
nChan = size(EEG.data,1);

X = double(EEG.data(:,1+lagAmount_samples:end))';
Y = double(EEG.data(:,1:end-lagAmount_samples))';

% Y = double(EEG.data(:,1+lagAmount_samples:end))';
% X = double(EEG.data(:,1:end-lagAmount_samples))';


[A B R U V] = canoncorr(X,Y);


% figure; stem(R);
% figure; stem(R.^2);
% 
% figure; imagesc(A);
% figure; imagesc(B);
% figure; imagesc(A-B);

% % numComps2Delete = 0;
% if numComps2Delete > nChan | numComps2Delete < 0
%     error('You specified an unreasonable number of components to delete');
% end
% badComps = length(R)+1-numComps2Delete:length(R); %R organized from large to small. Delete components with small correlation


badComps = find((R.^2) < Rsq_thres);
disp(['Removing ',num2str(length(badComps)),' bad comps with auto-CCA']);

%not plotting figure right now %CG 9/29/21
%figure; stem((R.^2));hold on; stem(badComps,R(badComps).^2,'r'); xlabel('CCA comp num');ylabel('Rsq');

%just quick notes below simplifying math to make sure we keep track of
%transposes properly. X_col means X array with column vectors. A_col means
%A matrix from canoncorr which used column vectors originally. X_row and
%A_row are the equivalent versions when flipping back to EEGLAB row vector
%format. 
% 
% U_col = X_col*A_col
% U_row = (U_col)' = (X_col*A_col)' = A_col' * X_col' = A_col' * X_row;
% U_row = A_col' * X_row;

% A_zeroed = A;
% A_zeroed(:,badComps) = deal(0);
% 
% cleanChData = A_zeroed*X;
% 
% EEG_clean.data = cleanChData';
nComp = size(A,2);
EEG_clean.icasphere = eye(nChan);%row
EEG_clean.icaweights = A'; %colum
EEG_clean.icawinv = pinv(EEG_clean.icaweights*EEG_clean.icasphere);
pop_editoptions('option_computeica',1);
EEG_clean = eeg_checkset(EEG_clean,'ica');
if ~isempty(badComps)
    EEG_clean = pop_subcomp( EEG_clean, badComps, 0);
end


end

