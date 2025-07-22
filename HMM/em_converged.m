function [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if
%   |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f is log lik.
% threshold defaults to 1e-4.

if nargin < 3
    threshold = 1e-4;
end

converged = 0;
decrease = 0;

if loglik - previous_loglik < -1e-3 % allow for a little imprecision
    fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
    decrease = 1;
end

% The following stopping criterion is from Numerical Recipes in C p423
delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end