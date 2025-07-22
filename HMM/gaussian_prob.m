function p=gaussian_prob(x, m, C, use_log)
% p=gaussian(x,m,C);
%
% Evaluate the multi-variate density with mean vector m and covariance
% matrix C for the input vector x.
% Vectorized version: Here X is a matrix of column vectors, and p is
% a vector of probabilities for each vector.

if nargin<4, use_log = 0; end

d=length(m);

if size(x,1)~=d
    x=x';
end
N=size(x,2);

m=m(:);
M=m*ones(1,N);
denom=(2*pi)^(d/2)*sqrt(abs(det(C)));
invC = inv(C);
mahal=sum(((x-M)'*invC).*(x-M)',2);   % Chris Bregler's trick
if use_log
    p = -0.5*mahal - log(denom);
else
    numer=exp(-0.5*mahal);
    p=numer/denom;
end