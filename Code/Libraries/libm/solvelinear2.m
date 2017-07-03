% Inputs:
%   N: number of points at which to compute solution
%   n: dimension of each equation system
%   func gives each residual in a separate column; rows for points
function x = solvelinear2(func,N,n,varargin)

  J = zeros(N,n*n);
  z = zeros(N,n);
  f0 = feval(func,z,varargin{:});
  for i=1:n
    x = z;
    x(:,i) = ones(N,1);
    J(:,(i-1)*n+1:i*n) = (feval(func,x,varargin{:})-f0);
  end;

  J = J';
  f0 = f0';

  x = zeros(n,N);
  for i=1:N
    jac = reshape(J(:,i),n,n);
    x(:,i) = -jac\f0(:,i);
  end;
  x = x';