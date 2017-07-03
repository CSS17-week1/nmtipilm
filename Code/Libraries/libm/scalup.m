% function y = scalup(x,xmin,xmax);
% Linearly rescales every element of the vector x 
% from [-1,1]  to  [xmin,xmax]
% where x,xmin,xmax can be vectors as long as they are all of the same 
% dimension
%
% November 9 1998

function xu = scalup(x,xmin,xmax);

[r c]   = size(x);
a       = (xmin+xmax)'/2;
b       = (xmax-xmin)'/2;
xu      = a*ones(1,r) + ( b*ones(1,r) ).*x';
xu      = xu';

% **********************************************************************

% **********************************************************************
