%function xd    = scaldown(x,xmin,xmax);
% Linearly scale a variable from [min,xmax] to [-1,1],
% where x,xmin,xmax can be vectors as long as they are all of the same 
% dimension
%
% November 9 1998
%

function xd     = scaldown(x,xmin,xmax);

[r c]   = size(x);
a       = 2*ones(r,c) ./ ( xmax - xmin );
b       = ones(r,c) - 2 * xmax ./ ( xmax - xmin );

xd      = b + a .* x;

% **********************************************************************

% **********************************************************************
