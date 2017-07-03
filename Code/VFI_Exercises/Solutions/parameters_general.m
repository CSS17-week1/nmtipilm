%%  set parameter values

num_realiz = 4;             % number of shock realizations
alpha  = 0.40;              % production parameter
beta   = 0.95;              % subjective discount factor

prob   = 1/num_realiz + zeros(num_realiz,num_realiz);
% prob(i,j) = prob(A(t+1)=Aj | A(t) = Ai)

delta  = .90;               % 1 - depreciation rate
A_high = 1.5;               % high value for technology
A_low  = 0.5;               % low value for technology
A_shock = linspace(A_low,A_high,num_realiz)';
% vector of technology shocks
convcrit = 1e-7;            % convergence criterion (epsilon)
rng('default');   % set the random generator seed
kmark = 10;  % set the initial value of capital for the simulation



%%   generate capital grid

mink =   0.01;                      % minimum value of the capital grid
maxk =  25.01;                      % maximum value of the capital grid
nk   = 1000;                        % number of grid points
kgrid = linspace(mink,maxk,nk)';     % the grid (linearly spaced)
ink = kgrid(2) - kgrid(1);          % increments


