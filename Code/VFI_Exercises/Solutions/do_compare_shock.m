%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Comparing solution of RBC model with or without irreversible investment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;



%% Generate solution and simulation for reversible investment case
% load parameters and grid
parameters;
% modify parameters here:
A_high = 1.25;               
A_low  = 0.75;              

% solve the model via VFI
vfi_AM;

% store results in few new variables
v_rev = v;
decis_rev = decis;
controls_rev = controls;


%% Generate solution and simulation for irreversible investment case
% load parameters and grid
parameters;
% modify parameters here:
A_high = 1.25;               
A_low  = 0.75;               

% solve the model via VFI
vfi_AM_irrinv;

% store results in few new variables
v_irrev = v;
decis_irrev = decis;
controls_irrev = controls;



%% FIGURES

figure(1)
plot(kgrid',v_rev(:,1),'r-',kgrid',v_rev(:,2),'r:',...
    kgrid',v_irrev(:,1),'b-',kgrid',v_irrev(:,2),'b:');
title('STOCH GROWTH MODEL: VALUE FUNCTION');
legend('rev_H', 'rev_L','irrev_H','irrev_L');
print value.ps

figure(2)
plot(kgrid',decis_rev(:,1),'r.',kgrid',decis_rev(:,2),'r:',kgrid',kgrid','k-',...
    kgrid',decis_irrev(:,1),'b.',kgrid',decis_irrev(:,2),'b:');
title('STOCH GROWTH MODEL: POLICY FUNCTION');
legend('rev_H', 'rev_L','45^o line','irrev_H','irrev_L');
axis([ 0 maxk 0 maxk ]);
print policy.ps

figure(3)
plot((1:n-1)',controls_rev(:,1),'r',...
    (1:n-1)',controls_irrev(:,1),'b');
title('STOCH GROWTH MODEL: CONSUMPTION');
legend('rev','irrev');
print consum.ps

figure(4)
plot((1:n-1)',controls_rev(:,2),'r',...
    (1:n-1)',controls_irrev(:,2),'b');
title('STOCH GROWTH MODEL: INVESTMENT');
legend('rev','irrev');
print invest.ps
