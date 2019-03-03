function y = dirichletrnd(alpha, n)
% DIRICHLETRND Sample from a Dirichlet distribution
%
% y = dirichletrnd(alpha)
% y = dirichletrnd(alpha, n)
%
% alpha: n-by-d matrix of parameters for Dirichlet-distribution. If alpha
%        is a 1-by-d matrix then n samples with the same parameters are
%        taken (Default: n=1).
%
% y: Samples from a Dirichlet distribution with parameter vector alpha.

% $Id: dirichletrnd.m 4824 2014-12-07 21:56:53Z johanl $

if nargin<2 || isempty(n), n=1; end

%replicated alpha
if size(alpha,1)==1 && n>1, alpha=repmat(alpha,[n 1]); end

%sample using gamrnd
y = gamrnd(alpha, 1);
y = bsxfun(@rdivide, y, sum(y,2));