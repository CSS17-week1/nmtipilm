% plot allocations as a function of the grid 
gridpoints = 1000; 
grid = linspace(phi_min,phi_max,gridpoints)';

g1_grid = s(1)*ones(gridpoints,1);
g2_grid = s(2)*ones(gridpoints,1);

a1_grid  = funeval(para1,fspace,   [grid   ]);
a2_grid  = funeval(para2,fspace,   [grid   ]);

lambda1_grid  = funeval(parlambda1,fspace,   [grid   ]);
lambda2_grid  = funeval(parlambda2,fspace,   [grid   ]);

prob1_grid_g1 =a1_grid.^nu;
prob1_grid_g2 =1 - prob1_grid_g1;

% with persistence
prob2_grid_g2 =perst_weight.*(a2_grid.^nu2) + (1-perst_weight).*(1-a2_grid.^nu2);
prob2_grid_g1 =1 - prob2_grid_g2;


dprob1_grid_g1 =nu.*(a1_grid.^(nu-1));
dprob1_grid_g2 = - dprob1_grid_g1;

% with persistence
dprob2_grid_g2 =perst_weight.*(nu2.*(a2_grid.^(nu2-1))) + (1-perst_weight).*(-nu2.*(a2_grid.^(nu2-1)) ) ;
dprob2_grid_g1 = - dprob2_grid_g2;


phiNext1_grid_g1_g1  =  grid + lambda1_grid.*(dprob1_grid_g1./prob1_grid_g1 );
phiNext1_grid_g2_g1  =  grid + lambda2_grid .*(dprob2_grid_g1./prob2_grid_g1 );
phiNext1_grid_g1_g2  =  grid + lambda1_grid.*(dprob1_grid_g2./prob1_grid_g2  );
phiNext1_grid_g2_g2  =  grid + lambda2_grid .*(dprob2_grid_g2./prob2_grid_g2  );

c1_grid  = ((grid).^(1/sig)); 


exp_planner_disc_utility_g1 = funeval(parexp_planner1,fspace ,  [grid  ]);
exp_planner_disc_utility_g2 = funeval(parexp_planner2,fspace ,  [grid   ]);


exp_disc_utility_g1 = funeval(parexp1,fspace ,  [grid   ]);
exp_disc_utility_g2 = funeval(parexp2,fspace ,  [grid   ]);

transfer_grid_g1 = c1_grid + g1_grid -1;
transfer_grid_g2 = c1_grid + g2_grid -1;
