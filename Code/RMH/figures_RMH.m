%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%     "Repeated Moral Hazard and Recursive Lagrangeans"   %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figures on the grid
plot_grid_results_RMH;


%%%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 1: Plot of the policy functions on the grid/1  %
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(2,2,1);
plot(grid,c1_grid,'b','LineWidth',2);
title('consumption'); xlabel(' \phi')
% legend('c1 low','c1 high');
subplot(2,2,2);
plot(grid,a1_grid,'b','LineWidth',2);
title('effort'); xlabel(' \phi')
% legend('a1 low','a1 high');
subplot(2,2,3);
plot(grid,phiNext1_grid_g1_g1,'b',grid,phiNext1_grid_g1_g2,'r','LineWidth',2);
hold on;
plot(grid,grid,'k:','LineWidth',.1);
title('\phi^H , \phi^L'); xlabel(' \phi')
hold off;
legend('\phi^H ','\phi^L','\phi');
subplot(2,2,4);
plot(grid,lambda1_grid,'b','LineWidth',2);
title('lambda'); xlabel(' \phi')
% legend('lambda1 low','lambda1 high');


%%%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 2: Plot of the policy functions on the grid/2  %
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(2,2,1);
plot(grid,transfer_grid_g1,'b',grid,transfer_grid_g2,'r','LineWidth',2);
title('transfers'); xlabel(' \phi');
legend('\tau^H','\tau^L');
subplot(2,2,2);
plot(grid,exp_planner_disc_utility_g1,'b',grid,exp_planner_disc_utility_g2,'r','LineWidth',2);
title('J(y,\phi)'); xlabel(' \phi')
legend('J(y^H,\phi)','J(y^L,\phi)');
subplot(2,2,3);
plot(grid,exp_disc_utility_g1,'b','LineWidth',2);
title('U(y,\phi)'); xlabel(' \phi')


% Simulate series
simul_RMH;

%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 3: simulated series, avg. allocations/1   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
subplot(2,2,1);
plot(1:periods_simulations,c_simu,'b','LineWidth',2);
title('average consumption'); xlabel('t')
subplot(2,2,2);
plot(1:periods_simulations,a_simu,'b','LineWidth',2);
title('average effort');xlabel('t')
subplot(2,2,3);
plot(1:periods_simulations,exp_disc_util_simu,'b','LineWidth',2);
title('U(y,\phi): average lifetime utility (agent)');xlabel('t')
subplot(2,2,4);
plot(1:periods_simulations,lambda_simu,'b','LineWidth',2);
title('average \lambda');xlabel('t')


%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 4: simulated series, avg. allocations/2   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
subplot(2,2,1);
plot(1:periods_simulations,phi1_simu(1:end-1),'b','LineWidth',2);
title('average \phi');xlabel('t')
subplot(2,2,2);
plot(1:periods_simulations,transfer_simu,'b','LineWidth',2);
title('average transfer');xlabel('t')
subplot(2,2,3);
plot(1:periods_simulations,exp_planner_disc_util_simu,'b','LineWidth',2);
title('J(y,\phi): average planner value');xlabel('t')
subplot(2,2,4);
plot(0:periods_simulations-2,bonds_simu(1:end-2),'b','LineWidth',2);%0:periods_simulations-2,bonds_simu2(1:end-2),'r',
title('average asset holdings');xlabel('t')


%%%%%%%%%%%%%%%%
%    FIGURE 5: the Pareto frontier %
%%%%%%%%%%%%%%%%

gam1_vec = linspace(phi_min,phi_max,1000);

figure(5);
plot(funeval(parexp1,fspace,gam1_vec'),funeval(parexp_planner1,fspace,gam1_vec') - ...
    gam1_vec'.*funeval(parexp1,fspace,gam1_vec'),'b-',...
    funeval(parexp1,fspace,gam1),funeval(parexp_planner1,fspace,gam1) - ...
    gam1.*funeval(parexp1,fspace,gam1),'ro','LineWidth',2);
title('Pareto Frontier');
xlabel('Agent expect. utility');
ylabel('Planner expect. utility');
