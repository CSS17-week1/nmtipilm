% function weightmc
% similar to weightx, but works in case
% where exogenous variables are described by finite Markov chain,
% with nExog states
%
%
% Inputs:
% z_current: K-vector of current exogenous states (with values in 1..nExog)
% TM: nExog x nExog Markov transition matrix
%   TM(i,j) gives probability of going from state j to state i
% X: (K times nExog)xD matrix of values over which to take expectations
%   the nExog consecutive numbers give the value of X in the nExog realizations
%   of next period;s shock
% Output:
% KxD matrix of expected values
% 
% written in C, compiled to DLL file
% Michael Reiter, Universitat Pompeu Fabra,  March 2004
% Feel free to use, copy and change
