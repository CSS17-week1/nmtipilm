%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model 
%   with collocation VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all; clc;


global alpha betta rho sig sigma delta sigeps ;
global nQuadr QuadrWeights QuadrPoints;
global  RoundAppr rounds_approx  fspace;


%% PARAMETERS:

alpha =.4;                     % production function coefficient
delta = .1;                    % depreciation rate for capital
betta =  .95;                   % discount factor
sig = 1; % .5;%                 % CRRA utility parameter (c^(1-sig))/(1-sig)
sigma = .05;                    % S.D. of productivity shock
rho = 0;%;.9; %                       % persistence of productivity shock
convcrit = 1e-8;                % convergence criterion

%% create a grid for  around det. SS
kstar = (1/(alpha*betta) - (1-delta)/alpha )^(1/(alpha-1));
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
rounds_approx = 2; % number of rounds of approximation
Order_vector = [ 5 10;  5 10]; % number of grid points for each round % 4 6 8
ntest = 100; % number of grid point for testing


%% Approximation type for CompEcon
%     approxtype = 'lin';  % piecewise linear
        approxtype = 'cheb'; % chebychev polynomials
%     approxtype = 'spli'; % splines
    splineorder = []; % splines' order, default are cubic splines

%% parameters for simulations
number_series =1;    %number of series we want generate           
periods_simulation =  100; % number of periods for the simulation
k0 = k_min.*ones(number_series,1);

%% Run main file 
mainVFISGM; % solves the model

% deliver an accuracy statistic (the smaller the better)
max_test

% save the solution
save solutionVFISGM.mat parvalue fspace ;

%% SIMULATIONS AND FIGURES

% simulate the series
[A, k,  c, y ,epsilon] = simul_VFISGM(number_series, ...
    periods_simulation, k0, parvalue);

figure(1);

subplot(2,1,1);
plot(1:periods_simulation, k);
xlabel('time'); ylabel('k'); %axis([1 periods_simulation 0 7]);
title('Simulation: capital'); legend('k_t');

subplot(2,1,2);
plot(1:periods_simulation-1,c(2:end));
xlabel('t'); ylabel('c'); 
title('Simulation: consumption'); legend('c_t');
