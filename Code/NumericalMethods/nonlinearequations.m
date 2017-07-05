myfun = inline('exp(-x.^2) - cos(x)');

%% bisection
x_bisect = bisect(myfun, -4, 6)
plot(-4:.01:6, myfun(-4:.01:6), 'b', ...
    -4:.01:6, zeros(length(-4:.01:6))); % Notice: at least 4 roots!

%% function iteration
myfung = inline('x - (exp(-x.^2) - cos(x))');
x_funiter = fixpoint(myfung, 0)

%% newton method

% create the function
fid = fopen ('newton_fun.m', 'w');
fprintf(fid, 'function [fval, fjac] = newton_fun(x);\n');
fprintf(fid, 'fval = exp(-x.^2) - cos(x);\n');
fprintf(fid, 'fjac = -2.*x.*exp(-x.^2) + sin(x) ;\n' );
fclose (fid);    

% perform the newton iterations
x_newton = newton('newton_fun', 4)

%% broyden's method
% solving the same function as all the other methods
x_broyden = broydn('newton_fun',4,1e-7,0,1)

% solving the intertemporal two-periods model
fid = fopen ('twoperiods.m', 'w');
fprintf(fid, 'function equ = twoperiods(c, s0);\n');
fprintf(fid, 'c1 = c(1); c2 = c(2);\n');
fprintf(fid, 'equ1 =  c1^-2 - .99 *(1.05)*(c2^-2);\n' );
fprintf(fid, 'equ2 =  c1 + c2/(1.05) - 1 - 1/(1.05) - s0;\n' );
fprintf(fid, 'equ =  [equ1; equ2];\n' );
fclose (fid);    

s0_vec = 0:0.1:10;

for i = 1:length(s0_vec)
    c_broyden(:, i) = broydn('twoperiods',[1; 1],1e-7,0,1, s0_vec(i) );
end
c_broyden




