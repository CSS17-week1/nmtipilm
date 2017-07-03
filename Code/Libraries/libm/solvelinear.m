% solves func(x)=0, assuming that func is linear in x
% Inputs:
%   func: handle or name of function, R^n->R^n
%   n: dimension of domain and range of func
%   varargin: will be passed to func
function x = solvelinear(func,n,varargin)

  J = zeros(n,n);
  z = zeros(n,1);
  f0 = feval(func,z,varargin{:});
  for i=1:n
    x = z;
    x(i) = 1;
    J(:,i) = feval(func,x,varargin{:})-f0;
  end;

  x = -J\f0;
