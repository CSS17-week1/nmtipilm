%% Exercise 3.1: Integrals

% number of nodes
n_vec = [11 21 31 101];

for i=1:length(n_vec)
    n = n_vec(i);

    % Calculate the integrals on the interval $[-1,1]$ using the trapezoid rule.
    [x,w] = qnwtrap(n, -1, 1);
    integral1(i) = w'*exp(-x);
    integral2(i) = w'*(abs(x).^0.5);
    integral3(i) = w'*(1./(1+ 25*x.^2));

    comparison_matrix_trap(i,:) = [n integral1(i) integral2(i) integral3(i)];
            
    % Calculate the integrals on the interval $[-1,1]$ using Simpson's rule.
    [x,w] = qnwsimp(n, -1, 1);
    integral1(i) = w'*exp(-x);
    integral2(i) = w'*(abs(x).^0.5);
    integral3(i) = w'*(1./(1+ 25*x.^2));

    comparison_matrix_simp(i,:) = [n integral1(i) integral2(i) integral3(i)];
            
    % Calculate the integrals on the interval $[-1,1]$ using Gauss-Legendre method.	
    [x,w] = qnwlege(n, -1, 1);	
    integral1(i) = w'*exp(-x);
    integral2(i) = w'*(abs(x).^0.5);
    integral3(i) = w'*(1./(1+ 25*x.^2));

    comparison_matrix_lege(i,:) = [n integral1(i) integral2(i) integral3(i)];

end
truevalue1 = -exp(-1)+ exp(1);
truevalue2 = 4/3;
truevalue3 = 1/5 * (atan(5) - atan(-5));
comparison_matrix = [0 truevalue1 truevalue2 truevalue3 ;  comparison_matrix_trap ; comparison_matrix_simp ;comparison_matrix_lege ]
	 