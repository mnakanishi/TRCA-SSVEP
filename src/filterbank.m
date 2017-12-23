function y = filterbank(eeg, fs, idx_fb)
% Filter bank design for decomposing EEG data into sub-band components [1].
%
% function y = filterbank(eeg, fs, idx_fb)
%
% Input:
%   eeg             : Input eeg data
%                     (# of channels, Data length [sample], # of trials)
%   fs              : Sampling rate
%   idx_fb          : Index of filters in filter bank analysis
%
% Output:
%   y               : Sub-band components decomposed by a filter bank.
%
% Reference:
%   [1] X. Chen, Y. Wang, S. Gao, T. -P. Jung and X. Gao,
%       "Filter bank canonical correlation analysis for implementing a 
%       high-speed SSVEP-based brain-computer interface",
%       J. Neural Eng., vol.12, 046008, 2015.
%
% Masaki Nakanishi, 22-Dec-2017
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

if nargin < 2
    error('stats:test_fbcca:LackOfInput', 'Not enough input arguments.'); 
end

if nargin < 3 || isempty(idx_fb)
    warning('stats:filterbank:MissingInput',...
        'Missing filter index. Default value (idx_fb = 1) will be used.'); 
    idx_fb = 1;
elseif idx_fb < 1 || 10 < idx_fb
    error('stats:filterbank:InvalidInput',...
        'The number of sub-bands must be 0 < idx_fb <= 10.'); 
end

[num_chans, ~, num_trials] = size(eeg);
fs=fs/2;

passband = [6, 14, 22, 30, 38, 46, 54, 62, 70, 78];
stopband = [4, 10, 16, 24, 32, 40, 48, 56, 64, 72];
Wp = [passband(idx_fb)/fs, 90/fs];
Ws = [stopband(idx_fb)/fs, 100/fs];
[N, Wn]=cheb1ord(Wp, Ws, 3, 40);
[B, A] = cheby1(N, 0.5, Wn);

y = zeros(size(eeg));
if num_trials == 1
    for ch_i = 1:1:num_chans
        y(ch_i, :) = filtfilt(B, A, eeg(ch_i, :));
    end % ch_i
else
    for trial_i = 1:1:num_trials
        for ch_i = 1:1:num_chans
            y(ch_i, :, trial_i) = filtfilt(B, A, eeg(ch_i, :, trial_i));
        end % trial_i
    end % ch_i
end % if num_trials == 1