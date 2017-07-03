% x=chebpol(ord,p);
% Purpose
% Create an ord-th order Chebyshev Polynomial 
% if p	is a vector then the polynomial will be created based 
%	on that vector
% November 9 1998
%
% ----------------------------------------------------------------------
function x=chebpol(ord,p);

	[r c]	= size(p);
 	x	= ones(r,ord+1);	
	x(:,2)	= p;
	if ord >= 2;  % corrected 9.2.2002, M.R.
	   for i	= 3:ord+1;
		x(:,i)	= 2.*p.*x(:,i-1)-x(:,i-2);
           end;
	end;

% **********************************************************************

% **********************************************************************
