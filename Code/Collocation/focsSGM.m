function equ = focsSGM(park, Grid,fspace);
% FOCSSGM contains the first order conditions for the basic RBC model
%
% (c) 2017 Antonio Mele 
 
% declare global variables
global alpha betta sig rho delta A_bar
global nQuadr QuadrWeights QuadrPoints

% rename the extrema of the functional space, for convenience
LowerBound = fspace.a;
UpperBound = fspace.b;


% rename grid columns for convenience
k = Grid(:,1);
A = Grid(:,2);

% evaluate policy functions
knext = funeval(park,fspace, Grid);
fofk = exp(A).*(k.^alpha) + (1-delta).*k;
c = max(fofk - knext, zeros(length(Grid),1));


n = length(k);
% generate nQuadr replications of the Grid, one for each realization of shock:
Grid_knext = kron(knext,ones(nQuadr,1));

% Expected value of next period A, corresponding to Grid:
ExpA =  rho*A;
% all realizations of next A:
GridANext = kron(ExpA,ones(nQuadr,1)) + ...
    kron(ones(n,1),QuadrPoints);
% truncate it to state space:
GridANext = min(max(GridANext,LowerBound(2)),UpperBound(2));

% grid in the next period, used for calculating expecations 
GridNext = [Grid_knext GridANext];

% calculate variables at t+1
knextnext = funeval(park,fspace, GridNext);
fofknext = exp(GridANext).*(Grid_knext.^alpha) + (1-delta).*Grid_knext;
 cnext = max(fofknext - knextnext, zeros(length(Grid_knext),1));
mucnext = muc(cnext);
mpknext = mpk(GridNext);

% calculate expectations with quadrature
exp_mucnext = (QuadrWeights'*reshape(mpknext.*mucnext,nQuadr,n))';

% Equation to be solved: Euler equation
equ = (muc(c) - betta.*exp_mucnext)./muc(c);

% avoid strange solutions
if (any(cnext<0))  || (any(c<0)) || (any(knext<0))  || (any(knextnext<0)) 
    equ(1) = 1e100;
end;
