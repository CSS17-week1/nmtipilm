% min1dnonconvex Computes local minimum of univariate function 
% on interval via Golden Search
% WORKS POINTWISE FOR FUNCTIONS R^n->R^n!!
% calling syntax and options follow GOLDEN of Miranda/Fackler,
% except for additional parameter ngrid, that defines the fineness
% of the grid on which we search
%
% USAGE
%   [x,fval] = goldsvec(f,a,b,ngrid,P1,P2,...);
% INPUTS
%   f         : name of function of form fval=f(x)
%   a,b       : left, right endpoints of interval
%   ngrid     : split (a,b) into ngrid intervals to do grid search
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : local maximum of f
%   fval    : function value estimate
%
% USER OPTIONS (SET WITH OPTSET)
%   tol     : convergence tolerance
%
% Michael Reiter, February 2004
function [x1,f1] = min1dnonconvex(f,a,b,ngrid,varargin)

  tol = optget('goldsvec','tol',sqrt(eps));

  np = length(a);
  res = zeros(ngrid+1,np);

  for(i=0:ngrid)
    lambda = i/ngrid;
    x = lambda*b + (1-lambda)*a;
    fx = feval(f,x,varargin{:});
    res(i+1,:) = fx';
  end;
  [m,indx] = min(res);
  
  indx_a = max(indx-1,1)';
  indx_b = min(indx+1,np-1)';
  lambda_a = (indx_a-1)/ngrid;
  lambda_b = (indx_b-1)/ngrid;
  
  aa = lambda_a.*b + (1-lambda_a).*a;
  bb = lambda_b.*b + (1-lambda_b).*a;

  [x1,f1] = goldsvec(f,aa,bb,varargin{:});
