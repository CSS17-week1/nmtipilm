%% Main collocation loop
%
% (c) 2017 Antonio Mele
tic;
for oo = 1: rounds_approx
    RoundAppr = oo;     % this tells you at which round of approximation we are
        
    % Range on which we approximate the solution:
    LowerBound = [k_min A_min ];
    UpperBound = [k_max A_max];
    
    % Approximation order
    Order = Order_vector(:,oo);
    
    % for reference, need this for guess after round 1
    if oo >= 2
        fspace_old = fspace; % save fspace from previous round, to calculate a guess
    end

    % display information about the round of approximation and the grid
    disp(' '); disp(' ' );
    disp(sprintf('RoundAppr %d, # gridpoints = [%d %d]',...
        RoundAppr, Order(1), Order(2)));
    disp(' '); 
    

    %% Generate the space of basis functions by using CompEcon toolbox for function approximation (see Miranda-Fackler book)
    % the following lines generate the basis function space
    % we can choose among Chebychev polynomials, splines of different orders and piecewise linear functions
    % uncomment the line that corresponds to the approximation you want to use
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
    end;
    
    % the following commands create efficient gridpoints (for example, when using Chebychev polynomials, they choose the zeroes)
    nodes = funnode(fspace);
    Grid = gridmake(nodes);
    
    % Set initial conditions
    if (RoundAppr == 1)
        knext = Grid(:,1);
    else    % if we are at second approx round, we use the solution of the first round as initial conditions on the new larger grid
        knext = funeval(park, fspace_old, Grid);
    end;    
    
    % generate basis functions Basis at Grid :
    Basis = funbas(fspace, Grid);    
    
    % set initial value for parameters of the approximation
    park =  Basis\knext;
    
    % solve FOCs with Broyden method for nonlinear equations
    [park,info] =  broydn('focsSGM',park,1e-8,0,1,Grid,fspace);
    disp(sprintf(' info = %d',info)); % if info=0, everything went fine, o/w the Broyden algorithm didn't converge
    disp(sprintf('    '));   
    
end;

% save the computation time
toc;
time_computation = toc;

% Calculate time needed for solution
time_hours = fix(toc/3600);
time_minutes = fix((toc - time_hours*3600)/60);
time_seconds = fix(toc -  time_hours*3600 - time_minutes*60);

disp(sprintf('    '));disp(sprintf('Time for solving the model was    '));
disp(sprintf('%d hours     %d  minutes %d  seconds',time_hours,time_minutes,time_seconds));

%% T E S T I N G    

% Testing the approximation: we calculate the value of the residuals at off-grid points. We therefore first generate a new grid with many points, then calculate the residuals on this grid and check the value of the norm and the max. If small, we can be relatively sure that the accuracy obtained is good

% Create a large grid with nphi points
k_nodes_test = linspace(LowerBound(1),UpperBound(1),ntest)';
A_nodes_test = linspace(LowerBound(2),UpperBound(2),ntest)';
GridTest = gridmake(k_nodes_test,A_nodes_test);

%calculate residuals on the new grid with optimized parameters
test_residuals =  focsSGM(park,GridTest,fspace);

% generate statistics to measure the accuracy
norm_test  = norm(test_residuals);
max_test   = max(abs(test_residuals));

testing_max = max_test(end);
testing_norm = norm_test(end);



