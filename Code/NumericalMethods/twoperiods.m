function equ = twoperiods(c, s0);
c1 = c(1); c2 = c(2);
equ1 =  c1^-2 - .99 *(1.05)*(c2^-2);
equ2 =  c1 + c2/(1.05) - 1 - 1/(1.05) - s0;
equ =  [equ1; equ2];
