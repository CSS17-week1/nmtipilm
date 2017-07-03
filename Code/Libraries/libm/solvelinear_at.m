% solves func(x)=0, assuming that func is linear in x
% same as solvelinear, but evaluates at z with step dz
function x = solvelinear_at(func,z,dz,varargin)

  n = length(z);
  J = zeros(n,n);
  f0 = feval(func,z,varargin{:});
  for i=1:n
    x = z;
    x(i) = z(i) + dz;
    J(:,i) = (feval(func,x,varargin{:})-f0)/dz;
  end;

  x = z-J\f0;
