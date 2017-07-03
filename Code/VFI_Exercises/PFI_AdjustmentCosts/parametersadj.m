%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Matlab code to solve a RBC model 
%   with investment adjustment costs by PFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  set parameter values


zeta = .25;                  % adjustment cost parameter
alpha  = 0.40;              % production parameter
beta   = 0.95;              % subjective discount factor 
prob   = [ .5 .5; .5 .5];   % prob(i,j) = probability (A(t+1)=Aj | A(t) = Ai)
delta  = .90;               % 1 - depreciation rate
A_high = 1.5;               % high value for technology
A_low  = 0.5;               % low value for technology
convcrit = 1e-7;            % convergence criterion (epsilon)
rng('default');             % set the random generator seed
kmark = 10;                 % set the initial value of capital for the simulation


%%   generate capital grid
   
mink =   0.01;                      % minimum value of the capital grid
maxk =  25.01;                      % maximum value of the capital grid   
nk   = 1000;                        % number of grid points
kgrid = linspace(mink,maxk,nk)';     % the grid (linearly spaced)
ink = kgrid(2) - kgrid(1);          % increments


