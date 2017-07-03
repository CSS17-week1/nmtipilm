%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model 
%   with irreversible investment by VFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  create the utility function matrices such that, for zero or negative
%  consumption, utility remains a large negative number so that
%  such values will never be chosen as utility maximizing


kapp = repmat(kgrid,1,nk);
kap = repmat(kgrid,1,nk)';

% substitute for the irreversible investment constraint
kapp(kapp<delta.*kap)= delta.*kap(kapp<delta.*kap);



cons = kron( kap.^alpha,A_shock')  + kron(delta.*kap - kapp, ones(1,num_realiz));
clear kapp kap

cons(cons<=0) = NaN;
util =  log(cons);
util(isnan(util)) = -inf;
clear cons


%%  initialize some variables

v       = zeros(nk,num_realiz);   % initial guess for value function:
% set to zero for simplicity
decis   = zeros(nk,num_realiz);   % initial value for policy function
metric  = 10;               % initial value for the convergence metric
iter = 0;
tme = cputime;
[rs,cs] = size(util);


%%  MAIN LOOP: iterate on Bellman's equation and get the decision rules and
%%  the value function at the optimum


while metric > convcrit;
    
    
    Ev = kron(v*prob',ones(1,nk));
    
    
    [tv,tdecis]=max(util +beta*Ev);
    
    tdecis=reshape(tdecis,num_realiz, nk)' ;
    tv=reshape(tv,num_realiz, nk)' ;
    
    metric=max(max(abs(tv-v)./tv));
    v=   tv;
    decis= tdecis;
    iter = iter+1;
    metric_vector(iter) = metric;
    disp(sprintf('iter = %g ; metric = %e', iter,metric));
end;
disp('  ');
disp(sprintf('computation time = %f', cputime-tme));

% transform the decision index in capital choice
decis=(decis-1)*ink + mink;



%%   simulate the model


kgrid = [ (0:ink:maxk)' ];  % capital grid
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state
states   = zeros(n-1,2);
controls = zeros(n-1,2);
[chain,state] = markov(prob,n,s0);
for i = 1:n-1;
    
    kprime = decis(kmark,chain(i));
    invest = kprime - delta*k;
    cons   = A_shock(chain(i))*k^(alpha) - invest;
    kmark = tdecis(kmark,chain(i));
    
    states(i,:) = [ k chain(i) ];
    controls(i,:) = [ cons invest ];
    k = kprime;
end;

