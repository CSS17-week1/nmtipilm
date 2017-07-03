% Program to check existence and uniqueness in linear RE models
% Input:
%   g0,g1,psi,pi,div: as defined in Sims' paper
% Output:
%   eu: stable solution exists iff eu(1)==1
%       solution is unique iff eu(2)==1
%   eigvals: generalized eigenvalues
% Michael Reiter, November2005
function [eu,eigvals] = checkeu(g0,g1,psi,pi,div)

  [Lambda Omega q z v]=qz(g0,g1);
  % Note that q'*Lambda*z' = g0 etc.,
  %  NOT  q*Lambda*z = g0 (as the Matlab help says!)


  realsmall = 1e-10; %much stricter than Sims
  dLambda = abs(diag(Lambda));
  dOmega = abs(diag(Omega));
  if(any(dLambda<realsmall & dOmega<realsmall))
    % taken from Sims:
    disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
    eu=[-2;-2];
    return;
  end

  [Lambda Omega q z]=qzdiv(div,Lambda,Omega,q,z);
  dLambda = abs(diag(Lambda));
  dOmega = abs(diag(Omega));
  dLambda = max(dLambda,realsmall); %to avoid dividing by 0;
  n = size(g0,1);
  n_unstable = sum(dLambda<=realsmall | dOmega./dLambda>div);
  
  iStable = 1:n-n_unstable;
  iUnstable = n-n_unstable+1:n;
  
  q1=q(iStable,:);
  q2=q(iUnstable,:);
  q2pi=q2*pi;
  q2psi=q2*psi;

  iExist = rank(q2pi) == rank([q2psi q2pi]);
  iUnique = rank(q*pi) == rank(q2pi);

  eu = [iExist;iUnique];
  warning off;
  ev = diag(Omega./Lambda);
  % sort according to absolute value:
  ev = sortrows([ev abs(ev)],2);
  eigvals = ev(:,1);
  warning on;



