
function equ = RMH_focs(par00, Gridphi,fspace);
 
global alpha betta sig phibar   epsil  nu nu2 gam1;
global phi_min phi_max s  ;
global   phi1  RoundAppr 
global  gam1_2 output_y perst_weight


LowerBound = fspace.a;
UpperBound = fspace.b;

% rename parameters
par0 = reshape(par00,length(par00)/8,8 );

parlambda1 =                    par0(:,1);
para1 =                         par0(:,2);
parexp_planner1 =               par0(:,3);
parexp_planner2 =               par0(:,4);
parexp1 =                       par0(:,5);
parexp2 =                       par0(:,6);
parlambda2 =                    par0(:,7);
para2 =                         par0(:,8);

% evaluate policy functions
a1 =                            funeval(para1,fspace, Gridphi);
a2 =                            funeval(para2,fspace, Gridphi);
lambda1 =                       funeval(parlambda1,fspace, Gridphi);
lambda2 =                       funeval(parlambda2,fspace, Gridphi);


%rename grid
phi1 = Gridphi(:,1);


% consumption is obtained directly from FOC consumption
c1 = ((  phi1).^(1/sig));

% probabilities for state of natures
prob1_s1 =  a1.^nu;
prob1_s2 = 1 - prob1_s1;

%derivative of probabilities
dprob1_s1 =  nu.*(a1.^(nu-1));
dprob1_s2 =  - dprob1_s1;

% allowing for persistence : works only if perst_weight is not zero

prob2_s2 = perst_weight.*(a2.^nu2) +(1-perst_weight).*(1 - a2.^nu2);
prob2_s1 = perst_weight.*(1 - a2.^nu2) +(1-perst_weight).*(a2.^nu2);

dprob2_s2 =  perst_weight.*(nu2.*(a2.^(nu2-1))) + (1-perst_weight).*(-nu2.*(a2.^(nu2-1)));
dprob2_s1 =  perst_weight.*(-nu2.*(a2.^(nu2-1))) + (1-perst_weight).*(nu2.*(a2.^(nu2-1)));


% likelihood ratios
like_ratio1_s1 =  dprob1_s1./prob1_s1;
like_ratio1_s2 = dprob1_s2./prob1_s2;

like_ratio2_s1 =  dprob2_s1./prob2_s1;
like_ratio2_s2 = dprob2_s2./prob2_s2;


% derivatives of likelihood ratios
dlike_ratio1_s1 =  -nu./(a1.^2);%
dlike_ratio1_s2 = (-nu*(nu-1)*(a1.^(nu-2)) - nu*(a1.^(2*nu-2)))./(prob1_s2.^2) ;

% allowing for persistence
dlike_ratio2_s2 =  perst_weight.*( -nu2./(a2.^2)) + ...
    (1- perst_weight).*((-nu2*(nu2-1)*(a2.^(nu2-2)) - nu2*(a2.^(2*nu2-2)))./(prob2_s2.^2)) ;%
dlike_ratio2_s1 = perst_weight.*((-nu2*(nu2-1)*(a2.^(nu2-2)) -...
    nu2*(a2.^(2*nu2-2)))./(prob2_s2.^2)) +(1- perst_weight).*( -nu2./(a2.^2))  ;


% Future costate variables
phiNext1_s1 = phi1 + lambda1.*like_ratio1_s1;

phiNext1_s2 = phi1 + lambda1.*like_ratio1_s2;

phiNext2_s1 = phi1 + lambda2.*like_ratio2_s1;

phiNext2_s2 = phi1 + lambda2.*like_ratio2_s2;



% future consumption 
c1Next_s1 = (( phiNext1_s1).^(1/sig)); %
c1Next_s2 = (( phiNext1_s2).^(1/sig)); %

c2Next_s1 = (( phiNext2_s1).^(1/sig)); %
c2Next_s2 = (( phiNext2_s2).^(1/sig)); %

% future effort
a1Next_s1 = funeval(para1,fspace,phiNext1_s1);
a1Next_s2 = funeval(para1,fspace,phiNext1_s2 );

a2Next_s1 = funeval(para2,fspace,phiNext2_s1);
a2Next_s2 = funeval(para2,fspace,phiNext2_s2 );


%instantaneous utility for agent next period
if sig==1
    utilityNext_s1 = log(c1Next_s1) - alpha*(a1Next_s1.^epsil);
    utilityNext_s2 = log(c1Next_s2) - alpha*(a1Next_s2.^epsil);
else
    utilityNext_s1 = (c1Next_s1.^(1-sig))./(1-sig) - alpha*(a1Next_s1.^epsil);
    utilityNext_s2 = (c1Next_s2.^(1-sig))./(1-sig) - alpha*(a1Next_s2.^epsil);
end;


if sig==1
    utilityNext2_s1 = log(c2Next_s1) - alpha*(a2Next_s1.^epsil);
    utilityNext2_s2 = log(c2Next_s2) - alpha*(a2Next_s2.^epsil);
else

    utilityNext2_s1 = (c2Next_s1.^(1-sig))./(1-sig) - alpha*(a2Next_s1.^epsil);
    utilityNext2_s2 = (c2Next_s2.^(1-sig))./(1-sig) - alpha*(a2Next_s2.^epsil);
end;

% this is function J : planner value function
exp_planner_disc_utility1 = funeval(parexp_planner1,fspace,Gridphi);
exp_planner_disc_utility1_next_s1 = funeval(parexp_planner1,fspace,phiNext1_s1 );
exp_planner_disc_utility1_next_s2 = funeval(parexp_planner2,fspace,phiNext1_s2 );

exp_planner_disc_utility2 = funeval(parexp_planner2,fspace,Gridphi);
exp_planner_disc_utility2_next_s1 = funeval(parexp_planner1,fspace,phiNext2_s1);
exp_planner_disc_utility2_next_s2 = funeval(parexp_planner2,fspace,phiNext2_s2 );



% This is U : agent continuation value
exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
exp_disc_utility1_next_s1 = funeval(parexp1,fspace,phiNext1_s1 );
exp_disc_utility1_next_s2 = funeval(parexp1,fspace,phiNext1_s2 );


exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);
exp_disc_utility2_next_s1 = funeval(parexp2,fspace,phiNext2_s1);
exp_disc_utility2_next_s2 = funeval(parexp2,fspace,phiNext2_s2 );


%%%%%%%%%%%%%%
%            %
%  EQUATIONS %
%            %
%%%%%%%%%%%%%%



% FOC lambda
equ1 = epsil*alpha.*(a1.^(epsil-1)) -   (betta.*dprob1_s1.*exp_disc_utility1_next_s1 ...
    +  betta.*dprob1_s2.*funeval(parexp2,fspace,phiNext1_s2 ));

equ1bis = epsil*alpha.*(a2.^(epsil-1)) -   (betta.*dprob2_s1.*funeval(parexp1,fspace,phiNext2_s1 ) ...
    +  betta.*dprob2_s2.*exp_disc_utility2_next_s2)  ;


% dynamic equat. for exp_disc_utility1
if sig==1
    equ2 = exp_disc_utility1 - log(c1) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob1_s1.*exp_disc_utility1_next_s1 ...
        -  betta.*prob1_s2.*funeval(parexp2,fspace,phiNext1_s2 ); 

    equ3 = exp_disc_utility2 - log(c1) ...
        + alpha.*(a2.^epsil) ...
        -  betta.*prob2_s1.*funeval(parexp1,fspace,phiNext2_s1 )  ...
        -  betta.*prob2_s2.*exp_disc_utility2_next_s2;

else

    equ2 = exp_disc_utility1 - (c1.^(1-sig))./(1-sig) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob1_s1.*exp_disc_utility1_next_s1 ...
        -  betta.*prob1_s2.*funeval(parexp2,fspace,phiNext1_s2 );


    equ3 = exp_disc_utility2 - (c1.^(1-sig))./(1-sig) ...
        + alpha.*(a2.^epsil) ...
        -  betta.*prob2_s1.*funeval(parexp1,fspace,phiNext2_s1 )  ...
        -  betta.*prob2_s2.*exp_disc_utility2_next_s2 ;
end;




% FOC a
equ4 =  -epsil*alpha.*(a1.^(epsil-1)).*phi1 - epsil*(epsil-1)*alpha.*(a1.^(epsil-2)).*lambda1 ...
    + ( betta.*lambda1.*(prob1_s1.*dlike_ratio1_s1.*utilityNext_s1 ...
    + prob1_s2.*dlike_ratio1_s2.*utilityNext_s2) ...
    +   betta.*dprob1_s1.*exp_planner_disc_utility1_next_s1 ...
    +   betta.*dprob1_s2.*funeval(parexp_planner2,fspace,phiNext1_s2))  ;

equ4bis =  -epsil*alpha.*(a2.^(epsil-1)).*phi1 - epsil*(epsil-1)*alpha.*(a2.^(epsil-2)).*lambda2 ...
    +( betta.*lambda2.*(prob2_s1.*dlike_ratio2_s1.*utilityNext2_s1 + ...
    prob2_s2.*dlike_ratio2_s2.*utilityNext2_s2) ...
    +   betta.*dprob2_s1.*funeval(parexp_planner1,fspace,phiNext2_s1) ...
    +   betta.*dprob2_s2.*exp_planner_disc_utility2_next_s2) ;



% dynamic equation for planner value
if sig==1
    equ5 = exp_planner_disc_utility1 - phi1.*(log(c1) ...
        - alpha.*(a1.^epsil)) ...
        -(output_y - c1 - s(1))  ...
        +  epsil*alpha.*(a1.^(epsil-1)).*lambda1 ...
        - (betta.*prob1_s1.*exp_planner_disc_utility1_next_s1 ...
        + betta.*prob1_s2.*funeval(parexp_planner2,fspace,phiNext1_s2)) ;

    equ6 = exp_planner_disc_utility2 - phi1.*(log(c1) ...
        - alpha.*(a2.^epsil)) ...
        -(output_y - c1 - s(2))  ...
        +  epsil*alpha.*(a2.^(epsil-1)).*lambda2...
        - (betta.*prob2_s1.*funeval(parexp_planner1,fspace,phiNext2_s1) ...
        + betta.*prob2_s2.*exp_planner_disc_utility2_next_s2)  ;

 
else
    equ5 = exp_planner_disc_utility1 - phi1.*((c1.^(1-sig))./(1-sig) ...
        - alpha.*(a1.^epsil)) ...
        -(output_y - c1 - s(1)) ...
        +  epsil*alpha.*(a1.^(epsil-1)).*lambda1 ...
        -(betta.*prob1_s1.*exp_planner_disc_utility1_next_s1 ...
        + betta.*prob1_s2.*funeval(parexp_planner2,fspace,phiNext1_s2)) ;


    equ6 = exp_planner_disc_utility2 - phi1.*((c1.^(1-sig))./(1-sig) ...
        - alpha.*(a2.^epsil)) ...
        -(output_y - c1 - s(2))  ...
        +  epsil*alpha.*(a2.^(epsil-1)).*lambda2...
        - (betta.*prob2_s1.*funeval(parexp_planner1,fspace,phiNext2_s1) ...
        + betta.*prob2_s2.*exp_planner_disc_utility2_next_s2)  ;


end;



% collect all equations in a vector
equ = [equ1; equ2; equ3;  equ4; equ5; equ6; equ1bis; equ4bis] ;

% avoid the solution is strange
if(any(a1<0)) || (any(c1<0))  ||(any(a1>1))   || (any( phiNext1_s1<0)) || (any( phiNext1_s2<0)) ||...
        (any(a2<0))   ||(any(a2>1))   || (any( phiNext2_s1<0)) || (any( phiNext2_s2<0))
    equ(1) = 1e100;
end;

