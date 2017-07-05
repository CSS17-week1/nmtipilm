%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model 
%   with collocation VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for oo = 1: rounds_approx
    RoundAppr = oo;     % this tells you at which round of approximation we are
    
    
    
    % Range on which we approximate the solution:
    LowerBound = [k_min A_min ];
    UpperBound = [k_max A_max];
    
    Order = Order_vector(:,oo);
    if oo >= 2
        fspace_old = fspace;
    end
    
    
    disp(' '); disp(' ' );
    disp(sprintf('RoundAppr %d, # gridpoints = [%d %d]',...
        RoundAppr, Order(1), Order(2)));
    disp(' ');
    
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
    
    
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
    end;
    
    % the following create gridpoints
    nodes = funnode(fspace);
    Grid = gridmake(nodes);
    
    % Set initial conditions
    if (RoundAppr == 1)
        
        value = zeros(length(Grid),1);
		
    else    % if we are at second approx round, we use the solution of the first round
        % as initial conditions on the new larger grid
        value = funeval(parvalue, fspace_old, Grid);
    end;
    
    
    % generate basis functions Basis at Grid :
    Basis = funbas(fspace,Grid);
    
    
    
    % set initial value for parameters of the approximation
    parvalue =  Basis\value;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                               %
    %               M A I N   L O O P               %
    %                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Solve Bellman equation with nonlinear solver
    [parvalue,info] =  broydn('residuals_VFISGM',parvalue,1e-8,0,1,Grid);
    disp(sprintf(' info = %d',info)); % if info=0, everything went fine, o/w the Broyden algorithm didn't converge
    disp(sprintf('    '));
    
    value = funeval(parvalue,fspace,Grid);
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
k_nodes_test = linspace(LowerBound(1),UpperBound(1),ntest)';%
A_nodes_test = linspace(LowerBound(2),UpperBound(2),ntest)';%

GridTest = gridmake(k_nodes_test,A_nodes_test);

%calculate residuals on the new grid with optimized parameters
test_residuals =  residuals_VFISGM(parvalue,GridTest);

% generatee statistics to measure the accuracy
norm_test  = norm(test_residuals);
max_test   = max(abs(test_residuals));



testing_max = max_test(end);
testing_norm = norm_test(end);

