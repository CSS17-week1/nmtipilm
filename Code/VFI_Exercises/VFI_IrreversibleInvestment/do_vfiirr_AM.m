%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model 
%   with irreversible investment by VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% do file for vfi_AM.m
clear all

% load parameters and grid
parametersirr;

% run the code
vfiirr_AM;

% generate figures
figuresirr;