function [mu,var] = discr_mom(x,prob)
  p = prob(:)';
  x = x(:);

  mu = p*x;
  var = p*(x-mu).^2;
  