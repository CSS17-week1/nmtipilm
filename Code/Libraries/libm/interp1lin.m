% 1-dimensional interpolation, same calling syntax as Matlab's interp1
% does no range check on x0!
% uses lookup from Miranda/Fackler
function z = interp1lin(x,y,x0)
  loc = lookup(x,x0,3);
  xl = x(loc);
  xr = x(loc+1);
  fac = (x0-xl)./(xr-xl);
  z = (1-fac).*y(loc) + fac.*y(loc+1);