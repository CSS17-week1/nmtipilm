%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a RBC model 
%   with irreversible investment by PFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do file for pfiirr_AM.m
clear all

% load parameters and grid
parametersirr;

% run the code
pfiirr_AM;

% generate figures
figuresirr;