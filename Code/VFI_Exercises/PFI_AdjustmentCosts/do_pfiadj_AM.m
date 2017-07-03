%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a RBC model 
%   with investment adjustment costs by PFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% do file for vfi_AM.m
clear all

% load parameters and grid
parametersadj;

% run the code
pfiadj_AM;

% generate figures
figuresadj;