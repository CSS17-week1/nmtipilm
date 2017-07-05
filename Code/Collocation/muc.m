function z = muc(c)
% MUC calculates the CRRA marginal utility of consumption
	global sig	
	z = c.^(-sig);