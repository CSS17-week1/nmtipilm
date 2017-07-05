%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% residuals VFI with collocation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = residuals_VFISGM(parvalue, Grid)



global alpha betta sig rho delta A_bar fspace
% global LowerBound UpperBound
global nQuadr QuadrWeights QuadrPoints

LowerBound = fspace.a;
UpperBound = fspace.b;

value = funeval(parvalue, fspace,Grid);

[~, rhsbellman] =  goldsvec('objectiveVFISGM',...
            zeros(length(Grid),1),...
            exp(Grid(:,2)).*(Grid(:,1).^alpha) + (1-delta).*Grid(:,1),... 
            parvalue,Grid);%

rhsbellman = - rhsbellman; % Remember: most optimizers are minimizers!!!
Res = value - rhsbellman;


