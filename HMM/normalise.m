function [M, c] = normalise(M)
% NORMALISE Make the entries of a (multidimensional) array sum to 1
% [M, c] = normalise(M)

c = sum(M(:));
% Set any zeros to one before dividing
d = c + (c==0);
M = M / d;

%if c==0
%  tiny = exp(-700);
%  M = M / (c+tiny);
%else
% M = M / (c);
%end