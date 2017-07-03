%% do file for vfi_AM.m
clear all

% load parameters and grid
parameters;
nk   = 100;                        % number of grid points
kgrid = linspace(mink,maxk,nk)';     % the grid (linearly spaced)
ink = kgrid(2) - kgrid(1);          % increments


% run the code
vfi_AM;

% generate figures
figures;