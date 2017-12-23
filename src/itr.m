function [ itr ] = itr(n, p, t)
% Calculate information transfer rate (ITR) for brain-computer interface 
% (BCI) [2]
% function [ itr ] = itr(n, p, t)
% 
% Input:
%   n   : # of targets
%   p   : Target identification accuracy (0 <= p <= 1) 
%   t   : Averaged time for a selection [s]
%
% Output:
%   itr : Information transfer rate [bits/min] 
%
% Reference:
%   [1] M. Cheng, X. Gao, S. Gao, and D. Xu,
%       "Design and Implementation of a Brain-Computer Interface With High 
%        Transfer Rates",
%       IEEE Trans. Biomed. Eng. 49, 1181-1186, 2002.
% 
% Masaki Nakanishi, 22-Dec-2017
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

if nargin < 3
    error('stats:itr:LackOfInput', 'Not enough input arguments.'); end

if p < 0 || 1 < p
    error('stats:itr:BadInputValue',...
        'Accuracy need to be between 0 and 1.');
elseif p < 1/n
    warning('stats:itr:BadInputValue',...
        'The ITR might be incorrect because the accuracy < chance level.');
    itr = 0;
elseif p == 1
    itr = log2(n)*60/t;
else
    itr = (log2(n) + p*log2(p) + (1-p)*log2((1-p)/(n-1)))*60/t;
end