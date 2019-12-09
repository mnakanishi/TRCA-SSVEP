function model = train_sscor(eeg, fs, num_fbs)
% Training stage of the sum of squared correlations (SSCOR)-based 
% steady-state visual evoked potentials (SSVEPs) detection [1].
%
% function model = train_sscor(eeg, fs, num_fbs)
%
% Input:
%   eeg         : Input eeg data 
%                 (# of targets, # of channels, Data length [sample])
%   fs          : Sampling rate
%   num_fbs     : # of sub-bands
%
% Output:
%   model       : Learning model for tesing phase of the ensemble 
%                 SSCOR-based method
%     - traindata   : Training data decomposed into sub-band components 
%                     by the filter bank analysis
%                     (# of targets, # of sub-bands, # of channels, 
%                      Data length [sample])
%     - W           : Weight coefficients for electrodes which can be 
%                     used as a spatial filter.
%     - num_fbs     : # of sub-bands
%     - fs          : Sampling rate
%     - num_targs   : # of targets
%
% See also:
%   test_sscor.m
%
% Reference:
%   [1] G. R. Kumar and M. R. Reddy,
%       "Designing a Sum of Squared Correlations Framework for Enhancing SSVEP
%        Based BCIs",
%       IEEE Trans. Neural Syst. Rehabil. Eng., vol. 27, pp. 2044-2050, 2019.
%
% Kuan-Jung Chiang and Masaki Nakanishi, 25-Nov-2019
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu


if nargin < 2
    error('stats:train_trca:LackOfInput', 'Not enough input arguments.'); 
end

if ~exist('num_fbs', 'var') || isempty(num_fbs), num_fbs = 3; end

[num_targs, num_chans, num_smpls, ~] = size(eeg);
trains = zeros(num_targs, num_fbs, num_chans, num_smpls);
W = zeros(num_fbs, num_targs, num_chans);
for targ_i = 1:1:num_targs
    eeg_tmp = squeeze(eeg(targ_i, :, :, :));
    for fb_i = 1:1:num_fbs
        eeg_tmp = filterbank(eeg_tmp, fs, fb_i);
        trains(targ_i,fb_i,:,:) = squeeze(mean(eeg_tmp, 3));
        w_tmp = sscor(eeg_tmp);
        W(fb_i, targ_i, :) = w_tmp;
    end % fb_i
end % targ_i
model = struct('trains', trains, 'W', W,...
    'num_fbs', num_fbs, 'fs', fs, 'num_targs', num_targs);


function W = sscor(eeg)
% Sum of Squared Correlation (SSCOR). This script was written based on
% the reference paper [1].
%
% function W = sscor(eeg)
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
%   [1] G. R. Kumar and M. R. Reddy,
%       "Designing a Sum of Squared Correlations Framework for Enhancing SSVEP
%        Based BCIs",
%       IEEE Trans. Neural Syst. Rehabil. Eng., vol. 27, pp. 2044-2050, 2019.
%
% Kuan-Jung Chiang and Masaki Nakanishi, 25-Nov-2019
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

[num_chans, ~, num_trials]  = size(eeg);
X_n = squeeze(mean(eeg,3));
X_n = bsxfun(@minus, X_n, mean(X_n,2));
K_1 = chol(X_n * X_n');

G_T_G = zeros(num_chans);
for trial_i = 1:1:num_trials
    x_i = squeeze(eeg(:,:,trial_i));
    x_i = bsxfun(@minus, x_i, mean(x_i,2));
    K_i = chol(x_i * x_i');
    C_1i = X_n * x_i';
    G_i = inv(K_1)' * C_1i / K_i;
    G_T_G = G_T_G + G_i' * G_i;
end % trial_i

[v_1,~] = eigs(G_T_G, 1);
W = K_1 \ v_1;
