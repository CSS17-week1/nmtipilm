% example of using function polyindx: 
%
% i = polyindx([3 1 4],[5 5 5],itensor)
%
% means that
% the i-th component of the parameter vector
% of a polynomial approximation of order 5
% [complete polynomials if itensor = 0,
% tensor product polynomials otherwise;
% order 3 means cubic, etc.]
% in 3 dimensions (x,y,z)
% gives the location of the element
% x^3 y z^4
% (of course, the elements are in reality Chebyshev polynomials,
% not monomials)
