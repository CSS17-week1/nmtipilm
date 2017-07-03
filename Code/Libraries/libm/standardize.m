% output is input, each column standardized
% to mean 0, variance 1
function xb = standardize(x)
  n = size(x,1);
  m = mean(x);
  x = x - repmat(m,n,1);
  s = sqrt(var(x));
  xb = x./repmat(s,n,1);