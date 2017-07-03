%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model by VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIGURES

figure(1)
plot(kgrid',v(:,1),'-',kgrid',v(:,2),':');
title('RBC MODEL: VALUE FUNCTION');
print value.ps

figure(2)
plot(kgrid',decis(:,1),'.',kgrid',decis(:,2),':',kgrid',kgrid','-');
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

