
%% FIGURES

figure(1)
plot(kgrid',v);
title('RBC MODEL: VALUE FUNCTION');
print value.ps

figure(2)
plot(kgrid',decis);
title('RBC MODEL: POLICY FUNCTION');
axis([ 0 maxk 0 maxk ]);
print policy.ps

figure(3)
plot((1:n-1)',controls(:,1));
title('RBC MODEL: CONSUMPTION');
print consum.ps

figure(4)
plot((1:n-1)',controls(:,2));
title('RBC MODEL: INVESTMENT');
print invest.ps

figure(5)
plot((1:iter)',log10(metric_vector));
title('RBC MODEL: CONVERGENCE CRITERION');
print convergence.ps