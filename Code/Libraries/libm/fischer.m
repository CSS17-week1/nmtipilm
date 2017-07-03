%Fischer's function to solve complementarity problems
%of the form
%  x>=a -> f(x)>=0
%  x<=b -> f(x)<=0
%cf. Miranda/Fackler Section 3.8
%we assume b>a.
%
%Inputs:
%f:        function value f(x)
%AminusX:  a-x
%BminusX:  b-x
%the function returns 0 if and only if
%1) f=0 & x>=a & x<=b
%OR
%2) f<=0 & x=a
%OR
%3) f>=0 & x=b
%
% corresponds to FOC of a MAXIMIZATION problem
function ftilde = fischer(f,AminusX,BminusX)
  aux = f + AminusX + sqrt(f.^2 + AminusX.^2);
  ftilde = aux + BminusX - sqrt(aux.^2 + BminusX.^2);
