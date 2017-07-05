%% Exercise 3.1: Integrals

% number of nodes
n = 11;

% Calculate the integral of $e^{-x}$ on the interval $[-1,1]$ using the trapezoid rule.
[x,w] = qnwtrap(n, -1, 1);
integral1 = w'*exp(-x);
	 	
% Calculate the integral of $\abs{x}^{\frac{1}{2}}$ on the interval $[-1,1]$ using Simpson's rule.
[x,w] = qnwsimp(n, -1, 1);
integral2 = w'*(abs(x).^0.5);
	 	
% Calculate the integral of $(1+25x^2)^{-1}$ on the interval $[-1,1]$ using Gauss-Legendre method.	
[x,w] = qnwlege(n, -1, 1);	
integral3 = w'*(1./(1+ 25*x.^2));

	 