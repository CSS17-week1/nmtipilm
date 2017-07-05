function [A, k,  c, y ,epsilon]  = simul_VFISGM(number_series,...
    periods_simulation, k0, parvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simul_VFISGM Simulates series for the stochastic growth model 
% solved with collocation on the Bellman equation
%
% USAGE: 
%   [A;k; c; y; epsilon] = simul(number_series, periods_simulation, k0, parvalue)
% INPUTS: 
%   number_series: number of series we want to create
%   periods_simulation: number of periods
%   k0: vector of initial conditions for capital
%   parvalue: interpolation coefficients
%
% simul_VFISGM returns simulated series for A
%   (productivity shock), k (capital), c (consumption), y (output), 
%   epsilon (the white noise shock in case A is autocorrelated). 
%   
%       Antonio Mele, October 2013, v. 0.1
%


global alpha delta rho sigeps

% allocate memory
A = zeros(number_series, periods_simulation);
k = zeros(number_series, periods_simulation);
c = zeros(number_series, periods_simulation);
y = zeros(number_series, periods_simulation);

% generate epsilon
epsilon = sigeps.*randn(number_series, periods_simulation);

% set initial conditions
A(:, 1) = epsilon(:, 1);
k(:,1) = k0;

% recursively generate series of interest
for ll=2:periods_simulation
    A(:,ll) = rho*A(:,ll-1) + epsilon(:,ll);
        
    k(:,ll) =  goldsvec('objectiveVFISGM',...
            zeros(number_series,1),...
            exp(A(:,ll)).*(k(:,ll-1).^alpha) + (1-delta).*k(:,ll-1),... 
            parvalue,[k(:,ll-1) A(:,ll)]);% 
    y(:,ll) = exp(A(:,ll)).*(k(:,ll-1).^(alpha));
    c(:,ll) = y(:,ll)+ (1- delta).*k(:,ll-1) - k(:,ll) ;
end

