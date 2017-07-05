function [A, k,  c, y ,epsilon]  = simul_SGM(number_series, periods_simulation, k0, park, fspace)

% simul_SGM Simulates series for the stochastic growth model solved with collocation
%
% USAGE: 
%   [A;k; c; y; epsilon] = simul_SGM(number_series, periods_simulation, k0, park, fspace)
% INPUTS: 
%   number_series: number of series we want to create
%   periods_simulation: number of periods
%   k0: vector of initial conditions for capital
%   park: interpolation coefficients
%   fspace: functional space (defined in the main file)
%
% simul_SGM returns simulated series for A (productivity shock), k (capital), c (consumption), y (output), epsilon (the white noise shock in case A is autocorrelated). 
%

% declare global variables
global alpha delta rho sigeps

% allocate memory
A = zeros(number_series, periods_simulation);
y = zeros(number_series, periods_simulation);
c = zeros(number_series, periods_simulation);
y = zeros(number_series, periods_simulation);

% generate epsilon
epsilon = sigeps.*randn(number_series, periods_simulation);

% set initial conditions
A(:, 1) = epsilon(:, 1);
k(:,1) = k0;

% recursively generate series of interest
for j=2:periods_simulation
    A(:,j) = rho*A(:,j-1) + epsilon(:,j);
    k(:,j) = funeval(park,fspace,[k(:,j-1) A(:,j)]); 
    y(:,j) = exp(A(:,j)).*(k(:,j-1).^(alpha));
    c(:,j) = y(:,j)+ (1- delta).*k(:,j-1) - k(:,j) ;
end
