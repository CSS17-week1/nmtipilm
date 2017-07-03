% function weightx
% multiplies a vector of weights with matrix
% helps to compute expected values efficiently
% 
% Inputs:
% W: Nx1 vector of (quadrature) weights
% X: (K times N)xD matrix
% 
% Output:
% KxD matrix, each element obtained by multiplying W with the
% corresponding subvector of X
% 
% written in C, compiled to DLL file
% Michael Reiter, Universitat Pompeu Fabra, 28.1.2003
% Feel free to use, copy and change
