% GJACOBI Solves Ax=b using Jacobi iteration
% USAGE
%   x = gjacobi(A,b,x);
% INPUTS
%   A     : nxn matrix
%   b     : n-vector
%   x     : n-vector of starting values, default=b
% OUTPUT
%   x     : approximate solution to Ax=b
%
% USER OPTIONS (SET WITH OPSET)
%   maxit : maximum number of iterations
%   tol   : convergence tolerence

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
% a bug fixed by Michael Reiter, 28.10.2005

function x = gjacobi(A,b,x)

if nargin<3, x=b; end
maxit = optget('gjacobi','maxit',1000);
tol   = optget('gjacobi','tol',sqrt(eps));

Q = diag(A);
for i=1:maxit
   dx = (b-A*x)./Q;  % error corrected Michael Reiter, 28.10.2005
   x = x + dx;
   if norm(dx)<tol, return; end
end
error('Maximum iterations exceeded in gjacobi')
