%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model by VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  create the utility function matrices such that, for zero or negative consumption, utility remains a large negative number so that such values will never be chosen as utility maximizing

cons1 = bsxfun(@minus, A_high*kgrid'.^alpha  + delta*kgrid' ,kgrid);
cons2 = bsxfun(@minus, A_low*kgrid'.^alpha  + delta*kgrid' ,kgrid);

cons1(cons1<=0) = NaN;
cons2(cons2<=0) = NaN;

util1 =  log(cons1);
util2 =  log(cons2);

util1(isnan(util1)) = -inf;
util2(isnan(util2)) = -inf;

%%  initialize some variables
v       = zeros(nk,2);   % initial guess for value function: set to zero for simplicity
decis   = zeros(nk,2);   % initial value for policy function
metric  = 10; % initial value for the convergence metric
iter = 0;
tme = cputime;
[rs,cs] = size(util1);

%%  MAIN LOOP: iterate on Bellman's equation and get the decision rules and the value function at the optimum
while metric > convcrit;    
    
    contv= beta*v*prob'; % continuation value
    
    [tv1,tdecis1]=max(bsxfun(@plus,util1,contv(:,1)) );  
    [tv2,tdecis2]=max(bsxfun(@plus,util2,contv(:,2)) );
    
    tdecis=[tdecis1' tdecis2'];
    tv=[tv1' tv2'];
    
    metric=max(max(abs((tv-v)./tv)));
    v=   tv; % .15*tv+.85*v; %
    decis= tdecis;%
    iter = iter+1;
    metric_vector(iter) = metric;
    disp(sprintf('iter = %g ; metric = %e', iter,metric));
end;
disp('  ');
disp(sprintf('computation time = %f', cputime-tme));

% transform the decision index in capital choice
decis=(decis-1)*ink + mink;

%%   simulate the model
kgrid = [ (0:ink:maxk)' ]; % capital grid
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state for the Markov chain
states   = zeros(n-1,2);
controls = zeros(n-1,2);
[chain,state] = markov(prob,n,s0); % generates the Markov chain
TFP = [A_high A_low]; % save productivity levels into TFP, for convenience

for i = 1:n-1
    kprime = decis(kmark, chain(i));
    invest = kprime - delta*k;
    cons = TFP(chain(i))*k^alpha - invest;
    kmark = tdecis(kmark, chain(i));
    states(i,:) = [ k chain(i) ];
    controls(i,:) = [ cons invest ];
    k = kprime;      
end

