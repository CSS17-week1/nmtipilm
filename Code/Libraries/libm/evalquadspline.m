% evaluate quadratic spline
% Inputs:
%   x,f,df: output from schumaker_convex
%   x0:     points at which_ to evaluate:
% Outputs:
%   f0,df0:  function_ values and derivatives at x0
% Michael Reiter, UPF, August 2006
function [f0,df0] = evalquadspline(x,f,df,x0)
  ilow = hunttable(x,x0);
  iupp = ilow+1;

  a = x(ilow);
  b = x(ilow+1);
  fa = f(ilow);
  dfa = df(ilow);
  dfb = df(ilow+1);
  


  d2fa = (dfb-dfa)./(b-a);
  d = x0-a;
  d2fad = d2fa.*d;
  df0 = dfa + d2fad;
  f0 = fa + d.*(dfa + 0.5*d2fad);

