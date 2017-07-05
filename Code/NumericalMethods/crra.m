function fval = crra(c1, s0);
c2 =(1 + s0 - c1)*1.05 + 1 ;
fval  = -(-(c1^-1 ) - .99 *(c2^-1 ));
