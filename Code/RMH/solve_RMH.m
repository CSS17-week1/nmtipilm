%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code solving and generating figures for the RMH solved with
% recursive Lagrangeans
%
%                       Antonio Mele
%                       March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;


global alpha betta sig phibar   epsil  nu nu2 gam1;
global phi_min phi_max s  ;
global   phi1  RoundAppr 
global  output_y perst_weight rounds_approx OrderVector Order_vector 

% PARAMETERS:

nba = 2;                        % number of states of the world
alpha =.50;                     % disutility of effort is alpha*v(a)
betta =  .95;                   % discount factor
sig = 2;                        % CRRA utility parameter (c^(1-sig))/(1-sig)
epsil = 2;                      % v(a) = a^epsil
output_y =1;                    % endowment in good times (in bad times is zero)
nphi = 100000;                     % number of gridpoints used for testing accuracy
perst_weight =0;                % parameter that governs the persistent of endowment: 
                                          %  if 0 no persistence, if 1 there is persistence

nu   =.5;                       % parameter of the prob(y_t+1= y^H| a_t) = a^nu;
nu2 = .5;                       % I allow for the possibility that transition probabilities are
                                % different in different states of the world

% outside option: The initial Pareto weight for the agent
gam1 =.5955;
% gam1_2 = [gam1; gam1];

% States of nature 
s = zeros(nba,1);
s(1) =  0;
s(2) =output_y;


% parameters for simulations
random_generations = 1;               % random_generations multiplied by metit gives
metit = 1;                                           % the total number of draws used; just for 
                                                           % technical reasons, we store
                                                           % the draws on big matrices
                                                           % with no more than 1000 rows.     
                                                           % For simulating one series,
                                                           % set both to 1
periods_simulations =  200; % number of periods for the simulation


rounds_approx = 2; % number of rounds of approximation
Order_vector = [10 30]; % number of grid points for each round
mainRMH; % solves the model

figures_RMH; % generates figures

