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
