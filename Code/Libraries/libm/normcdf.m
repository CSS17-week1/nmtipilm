function p = normcdf(x,mu,sigma)
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end
p = 0.5 * erfc(-(x-mu)./(sqrt(2)*sigma));
