% GOLDSVEC Computes local minimum of univariate function 
% on interval via Golden Search
% WORKS POINTWISE FOR FUNCTIONS R^n->R^n!!
%
% USAGE
%   [x,fval] = goldsvec(f,a,b,P1,P2,...);
% INPUTS
%   f         : name of function of form fval=f(x)
%   a,b       : left, right endpoints of interval
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : local maximum of f
%   fval    : function value estimate
%
% USER OPTIONS (SET WITH OPTSET)
%   tol     : convergence tolerance

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x1,f1] = goldsvec(f,a,b,varargin)

tol = optget('goldsvec','tol',sqrt(eps));

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
d  = b-a;
x1 = a+alpha1*d;
x2 = a+alpha2*d;
f1 = feval(f,x1,varargin{:});
f2 = feval(f,x2,varargin{:});

xx = 0*x1;

d = alpha1*alpha2*d;
while any(d>tol)
  d = d*alpha2;
  if2bigger = f2>f1;
  i1 = find(if2bigger);
  i2 = find(if2bigger-1);
  x2(i1) = x1(i1); 
  x1(i1) = x1(i1)-d(i1); 
  f2(i1) = f1(i1); 
  x1(i2) = x2(i2); 
  x2(i2) = x2(i2)+d(i2); 
  f1(i2) = f2(i2); 
  
  xx(i1) = x1(i1);
  xx(i2) = x2(i2);
  ff = feval(f,xx,varargin{:});
  f1(i1) = ff(i1);
  f2(i2) = ff(i2);
end

% Return the smaller of the two

i2less = find(f2<f1);
x1(i2less) = x2(i2less);
f1(i2less) = f2(i2less); 
