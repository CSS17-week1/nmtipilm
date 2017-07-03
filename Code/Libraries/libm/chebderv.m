% gives the coefficient of the derivative of a Chebyshev polynomial
% Inputs:
%    c: coefficien of the Chebyshev polynomial, c(1) is the constant, etc.
%    a: lower bound 
%    b: lower bound 
% Output:
%    Cheb. coefficients of the derivative, same size as input c
%       (means that last element is set to zero)
function cder = chebderv(c,a,b)
  n = size(c,1);
  cder = zeros(n,1);
  cder(n-1)=2*(n-1)*c(n);
  for (j=n-3 : -1 : 0)
    cder(j+1)=cder(j+3)+2*(j+1)*c(j+2);
  end;
  con=2.0/(b-a);
  for (j=0:n-1)
    cder(j+1) = cder(j+1) * con;
  end;
  cder(1) = cder(1) / 2;  

