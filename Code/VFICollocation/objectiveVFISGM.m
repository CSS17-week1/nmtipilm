function z = objectiveVFISGM(kprime, parvalue, Grid)

global alpha betta sig rho delta A_bar fspace
% global LowerBound UpperBound
global nQuadr QuadrWeights QuadrPoints

LowerBound = fspace.a;
UpperBound = fspace.b;


%rename grid
k = Grid(:,1);
A = Grid(:,2);

if sig==1
	u = log(exp(A).*(k.^alpha) +(1-delta).*k - kprime);
else
	u = ((exp(A).*(k.^alpha) +(1-delta).*k - kprime).^(1-sig))./(1-sig);
end

n = size(Grid,1);
% generate nQuadr replications of the Grid, one for each realization of shock:
Grid_kprime = kron(kprime,ones(nQuadr,1));

% Expected value of next period A, corresponding to Grid:
ExpA =  rho*A;
% all realizations of next A:
GridANext = kron(ExpA,ones(nQuadr,1)) + ...
    kron(ones(n,1),QuadrPoints);
% truncate it to state space:
GridANext = min(max(GridANext,LowerBound(2)),UpperBound(2));
GridNext = [Grid_kprime GridANext];
% calculate the value function at t+1
valuenext = funeval(parvalue, fspace,GridNext); 
% calculate the expectations
expectedvalue = (QuadrWeights'*reshape(valuenext,nQuadr,n))';

z1 = u + betta.*expectedvalue; 
z = -z1;
