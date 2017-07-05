function [fval, fjac] = newton_fun(x);
fval = exp(-x.^2) - cos(x);
fjac = -2.*x.*exp(-x.^2) + sin(x) ;
