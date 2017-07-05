%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Matlab code for solving stochastic growth model 
%  with orthogonal collocation over FOCs (DO FILE)
%
% (c) 2017 Antonio Mele
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
close all; clear all; clc;

% Declare global variables
global alpha betta rho sig sigma delta sigeps ;
global nQuadr QuadrWeights QuadrPoints;
global  RoundAppr rounds_approx  ;

%% PARAMETERS:

alpha =.4;                     % production function coefficient
delta = .1;                    % depreciation rate for capital
betta =  .95;                   % discount factor
sig = 1; % .5;%              % CRRA utility parameter (c^(1-sig))/(1-sig)
sigma = .05;                    % S.D. of productivity shock
rho = 0;%;  .9; %                     % persistence of productivity shock

%% create a grid for capital
kstar = (1/(alpha*betta) - ...
    (1-delta)/alpha )^(1/(alpha-1)); % deterministic steady state for capital
k_min = .5*kstar;
k_max = 2*kstar;

%% parameters for quadrature 
nQuadr = 50; % 100;% %number of quadrature points;
% we choose nQuadr high to get smoothness;
[QuadrPoints,QuadrWeights] = qnwnorm(nQuadr,0,sigma^2);

%% Range for shock
sigeps = sigma/sqrt(1-rho^2);
% Range for shock:
A_max = 3*sigeps;
A_min =  -3*sigeps;


%% Parameters for the collocation algorithm
rounds_approx = 2; % number of approximation rounds
Order_vector = [ 5 10;  5 10]; % number of grid points for each round 
ntest = 100; % number of grid points for testing (for each dimension)


%% Approximation type for CompEcon
% approxtype = 'lin';  % piecewise linear
approxtype = 'cheb'; % chebychev polynomials
% approxtype = 'spli'; % splines
splineorder = []; % splines' order, default are cubic splines

%% parameters for simulations
number_series = 1;    % number of series          
periods_simulation =  100; % number of periods for the simulation
k0 = k_min.*ones(number_series,1); % initial value for capital in the simulations

%% Run main file 
mainSGM; % solves the model

% deliver an accuracy statistic (the smaller the better)
max_test

% save the solution
save solutionSGM.mat park fspace ;

%% SIMULATIONS AND FIGURES
% simulate the series
[A, k,  c, y ,epsilon] = simul_SGM(number_series, ...
    periods_simulation, k0, park, fspace);

figure(1); % capital
subplot(2,1,1);
plot(1:periods_simulation, k);
xlabel('time'); ylabel('k'); axis([1 periods_simulation 0 7]);
title('Simulation: capital'); legend('k_t');

subplot(2,1,2); % consumption
plot(1:periods_simulation-1,c(2:end));
xlabel('t'); ylabel('c'); 
title('Simulation: consumption'); legend('c_t');
