%% Golden search

% function of one variable
[x_golden, fval_golden] = golden('humps', -10, 10)
[x_golden2, fval_golden2] = golden('humps', .2, 2) % there are several local maximizers!

% function with several variables
% REMINDER: this is a minimizer!!!
fid = fopen ('crra.m', 'w');
fprintf(fid, 'function fval = crra(c1, s0);\n');
fprintf(fid, 'c2 =(1 + s0 - c1)*1.05 + 1 ;\n');
fprintf(fid, 'fval  = -(-(c1^-1 ) - .99 *(c2^-1 ));\n' );
fclose (fid);    

s0_vec = 0:0.1:10;
for i = 1:length(s0_vec)
    [x_goldsvec(i), fval_goldsvec(i)] = goldsvec('crra', 0, 1+ 1/1.05 + s0_vec(i), s0_vec(i));
end

[x_goldsvec' fval_goldsvec' s0_vec']