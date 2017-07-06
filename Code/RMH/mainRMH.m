%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
% This code solve a dynamic agency problem          %
% with the recursive Lagrangean approach            %
% and collocation method oveer Lagrangean FOCs      %
%                                                   %
%                                                   %
%                    Antonio Mele, March 2010     %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
for oo = 1: rounds_approx
    RoundAppr = oo;     % this tells you at which round of approximation we are

    % create a grid for costates around the first pareto weight for agent
    %     gam1 =    mean(gam1_2);%gam1_2(1,1);%
    phibar = gam1;
    phi_sd =.3*gam1;
    phi_min = phibar- phi_sd;
    phi_max = phibar+phi_sd;

    % Range on which we approximate the solution:
    LowerBound = [phi_min  ];
    UpperBound = [phi_max ];

    Order = Order_vector(oo);
    if oo >= 2
        fspace_old = fspace;
    end



    disp(sprintf('  RoundAppr = %d ',RoundAppr));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % generate the space of basis functions
    % by using CompEcon toolbox for function
    % approximation (see Miranda-Fackler book)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % the following lines generate basis function space
    % we can choose among chebychev polynomials, splines of different
    % orders and piecewise linear functions

    % approxtype = 'lin';  % piecewise linear
    % approxtype = 'cheb'; % chebychev polynomials
    approxtype = 'spli'; % splines
    splineorder = []; % splines' order, default are cubic splines
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
    end;

    % the following create gridpoints
    nodes = funnode(fspace);
    Gridphi = gridmake(nodes);

    % Load initial conditions from an external file. These condtions were
    % found by setting nu = 0 and solving the model iteratively while
    % increasing the value of nu till around .5
    if (RoundAppr == 1)

        load initial_conditions_RMH a1 lambda1 exp_planner_disc_utility1  exp_planner_disc_utility2 ...
            exp_disc_utility1  exp_disc_utility2 a2 lambda2

    else    % if we are at second approx round, we use the solution of the first round
        % as initial conditions on the new larger grid
        a1 =funeval(para1,fspace_old, Gridphi);
        a2 = funeval(para2,fspace_old, Gridphi)  ;
        lambda1 =funeval(parlambda1 , fspace_old, Gridphi)  ;
        lambda2 =funeval(parlambda2 , fspace_old, Gridphi)  ;
        exp_disc_utility1 = funeval(parexp1,fspace_old,Gridphi);
        exp_disc_utility2 = funeval(parexp2,fspace_old,Gridphi);
        exp_planner_disc_utility1  = funeval(parexp_planner1,fspace_old,Gridphi);
        exp_planner_disc_utility2  = funeval(parexp_planner2,fspace_old,Gridphi);

    end;


    % generate basis functions Chebbasis at Gridphi :
    Basis = funbas(fspace,Gridphi);


    % set the name of costate variable to phi1
    phi1 = Gridphi(:,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set initial value for parameters of the approximation
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    parlambda1 =  Basis\lambda1;
    parlambda2 =  Basis\lambda2;

    para1 =  Basis\a1;
    para2 =  Basis\a2;

    parexp_planner1 =  Basis\exp_planner_disc_utility1;
    parexp_planner2 =  Basis\exp_planner_disc_utility2;

    parexp1 =  Basis\exp_disc_utility1 ;
    parexp2 =  Basis\exp_disc_utility2 ;

    % save all parameters in the vector parpolicy
    parpolicy = [parlambda1; para1 ; parexp_planner1 ; parexp_planner2 ;...
        parexp1 ; parexp2; parlambda2; para2];%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                               %
    %               M A I N   L O O P               %
    %                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solve Lagrangean FOCs by Broyden method for nonlinear equations
    [opt_vec,info] =  broydn('RMH_focs',parpolicy,1e-8,0,1,Gridphi,fspace);
    disp(sprintf(' info = %d',info)); % if info=0, everything went fine, o/w the Broyden algorithm didn't converge
    disp(sprintf('    '));

    %rename optimized parameters
    par00 = opt_vec;


    par0 = reshape(par00,length(Gridphi),8 );

    parlambda1 =                par0(:,1);
    para1 =                     par0(:,2);
    parexp_planner1 =           par0(:,3);
    parexp_planner2 =           par0(:,4);
    parexp1 =                   par0(:,5);
    parexp2 =                   par0(:,6);
    parlambda2 =                par0(:,7);
    para2 =                     par0(:,8);




    % calculate approximated policy functions on the grid
    a1 =                     funeval(para1,fspace, Gridphi);
    a2 =                     funeval(para2,fspace, Gridphi);

    lambda1 = funeval(parlambda1,fspace, Gridphi);
    lambda2 = funeval(parlambda2,fspace, Gridphi);

    exp_planner_disc_utility1 = funeval(parexp_planner1,fspace,Gridphi);
    exp_planner_disc_utility2 = funeval(parexp_planner2,fspace,Gridphi);

    exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
    exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);

    c1 = (( phi1).^(1/sig));



    % update vector of coefficients
    parpolicy = opt_vec;


end;


toc;
time_computation = toc;

% Calculate time needed for solution
time_hours = fix(toc/3600);
time_minutes = fix((toc - time_hours*3600)/60);
time_seconds = fix(toc -  time_hours*3600 - time_minutes*60);

disp(sprintf('    '));disp(sprintf('Time for solving the model was    '));
disp(sprintf('%d hours     %d  minutes %d  seconds',time_hours,time_minutes,time_seconds));





%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
%     T E S T I N G    %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%


% Testing the approximation: Residuals at points off the grid:

% Create a large grid with nphi points
grid_nodes_test = linspace(LowerBound,UpperBound,nphi)';%
GridTest = gridmake(grid_nodes_test);

%calculate residuals on the new grid with optimized parameters
test_residuals =  RMH_focs(parpolicy,GridTest,fspace);

% generatee statistics to measure the accuracy
norm_test  = norm(test_residuals);
max_test   = max(abs(test_residuals));



testing_max = max_test(end);
testing_norm = norm_test(end);



