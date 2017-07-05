function z = mpk(Grid)
% MPK calculates the marginal product of capital on Grid

	global alpha delta
	k = Grid(:,1);
    A = Grid(:,2);
	z = exp(A).*alpha.*(k.^(alpha-1)) + 1 - delta;
	