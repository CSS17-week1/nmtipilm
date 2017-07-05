% Univariate normal distribution

[x, w] = qnwnorm(30, 0,1);

sigmavec = [.5 1 2 10];

for i=1:length(sigmavec)

    expectations_uni(i) = w'*(x.^(-sigmavec(i)));

end
expectations_uni


% Multivariate normal distribution

[x, w] = qnwnorm([10, 15], [ 3 ; 4 ], [2 -1; -1 4 ]);

expectations_multi = w'*exp(x(:, 1) + x(:, 2))