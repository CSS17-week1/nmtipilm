%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Comparative static
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Generate solution and simulation for case A
% load parameters and grid
parameters;

% solve the model via VFI 
vfi_AM;

% store results in few new variables
v_A = v;
decis_A = decis;
controls_A = controls;


%% Generate solution and simulation for case B
% load parameters and grid
parameters;

% modify parameters here
prob   = [ .95 .05; .05 .95];   


% solve the model via VFI
vfi_AM;

% store results in few new variables
v_B = v;
decis_B = decis;
controls_B = controls;



%% FIGURES

figure(1)
plot(kgrid',v_A(:,1),'r-',kgrid',v_A(:,2),'r:',...
    kgrid',v_B(:,1),'b-',kgrid',v_B(:,2),'b:');
title('STOCH GROWTH MODEL: VALUE FUNCTION');
legend('A_H', 'A_L','B_H','B_L');
print value.ps

figure(2)
plot(kgrid',decis_A(:,1),'r.',kgrid',decis_A(:,2),'r:',kgrid',kgrid','k-',...
    kgrid',decis_B(:,1),'b.',kgrid',decis_B(:,2),'b:');
title('STOCH GROWTH MODEL: POLICY FUNCTION');
legend('A_H', 'A_L','45^o line','B_H','B_L');
axis([ 0 maxk 0 maxk ]);
print policy.ps

figure(3)
plot((1:n-1)',controls_A(:,1),'r',...
    (1:n-1)',controls_B(:,1),'b');
title('STOCH GROWTH MODEL: CONSUMPTION');
legend('A','B');
print consum.ps

figure(4)
plot((1:n-1)',controls_A(:,2),'r',...
    (1:n-1)',controls_B(:,2),'b');
title('STOCH GROWTH MODEL: INVESTMENT');
legend('A','B');
print invest.ps
