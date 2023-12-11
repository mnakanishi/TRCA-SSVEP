function W = ftrca(eeg)
% Fast task-related component analysis (fTRCA). This script was written based on
% the reference paper [1].
%
% function W = ftrca(eeg)
%
% Input:
%   eeg         : Input eeg data 
%                 (# of channels, Data length [sample], # of trials)
%
% Output:
%   W           : Weight coefficients for electrodes which can be used as 
%                 a spatial filter.
%   
% Reference:
%   [1] H. Tanaka, T. Katura, H. Sato,
%       "Task-related component analysis for functional neuroimaging and 
%        application to near-infrared spectroscopy data",
%       NeuroImage, vol. 64, pp. 308-327, 2013.
%
%   [2] K. -J. Chiang, C. M. Wong, F. Wan, T. -P. Jung, M. Nakanishi,
%       "Reformulating task-related component analysis to reduce its
%       computational cost",
%       Under review. 
%
% Chi Man Wong, Kuan-Jung Chiang, and Masaki Nakanishi, 08-Sep-2022 
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

[num_chans, num_smpls, num_trials]  = size(eeg);
for trial_i = 1:1:num_trials
    x1 = squeeze(eeg(:,:,trial_i));
    eeg(:,:,trial_i) = bsxfun(@minus, x1, mean(x1,2));
end % trial_i
SX = sum(eeg,3);
S = SX*SX';
UX = reshape(eeg, num_chans, num_smpls*num_trials);
Q = UX*UX';
[W,~] = eigs(S, Q);
