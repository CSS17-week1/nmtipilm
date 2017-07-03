% bisectvec: bisection method for one-dmensional root finding
% USAGE
%   x=bisect(f,a,b,P1,P2,...)
% INPUTS
%   func: name of function with
%          input:  array of any dimension
%          output: array of the same dimension
%   a: array of lower bounds
%   b: array of upper bounds
%   P1,P2,...: optional extra parameters to be passed to func
% OUTPUT
%   x: the roots of func, i.e., func(x)=0
% Michael Reiter, December 2005
function x=bisectvec(func,a,b,varargin)
  tol  = 1e-8;
  fa = feval(func,a,varargin{:});
  fb = feval(func,b,varargin{:});
  if any(fa.*fb>0)
    error('root not bracketed') 
  end
  while any(abs(b-a)>tol)
    x = (a+b)/2;
    f = feval(func,x,varargin{:});
    ib = f.*fa <=0;
    ia = ~ib;
    b(ib) = x(ib);
    fb(ib) = f(ib);
    a(ia) = x(ia);
    fa(ia) = f(ia);
  end
