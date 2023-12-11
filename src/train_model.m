function model = train_model(eeg, fs, num_fbs, algorithm)
% Training a model based on the template-matching method for 
% steady-state visual evoked potentials (SSVEPs) detection [1].
%
% function model = train_model(eeg, fs, num_fbs)
%
% Input:
%   eeg         : Input eeg data 
%                 (# of targets, # of channels, Data length [sample])
%   fs          : Sampling rate
%   num_fbs     : # of sub-bands
%
% Output:
%   model       : Learning model for tesing phase of the ensemble method
%     - traindata   : Training data decomposed into sub-band components 
%                     by the filter bank analysis
%                     (# of targets, # of sub-bands, # of channels, 
%                      Data length [sample])
%     - W           : Weight coefficients for electrodes which can be 
%                     used as a spatial filter.
%     - num_fbs     : # of sub-bands
%     - fs          : Sampling rate
%     - num_targs   : # of targets
%     - algorithm   : Spatial filtering method (ftrca, trca, or sscor)
%
% See also:
%   test_model.m
%
% Reference:
%   [1] M. Nakanishi, Y. Wang, X. Chen, Y. -T. Wang, X. Gao, and T.-P. Jung,
%       "Enhancing detection of SSVEPs for a high-speed brain speller using 
%        task-related component analysis",
%       IEEE Trans. Biomed. Eng, 65(1): 104-112, 2018.
%
% Masaki Nakanishi, 08-Sep-2022
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

if nargin < 4 || ~exist('algorithm', 'var')
  warning('stats:train_model:LackOfInput', 'Algorithm not specified. Deault model (fTRCA) will be used.');
  algorithm = 'ftrca';
elseif nonzeros(strcmp(algorithm, {'ftrca', 'trca', 'sscor'}))
  disp(['The algorithm, ' algorithm ', will be used.']);
else
  error('stats:train_model:WrongInput', 'An algorithm should be selected from ftrca, trca, or sscor.');
end

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
    eval(['w_tmp = ' algorithm '(eeg_tmp);'])
    W(fb_i, targ_i, :) = w_tmp(:,1);
  end % fb_i
end % targ_i
model = struct('trains', trains, 'W', W,...
    'num_fbs', num_fbs, 'fs', fs, 'num_targs', num_targs);

