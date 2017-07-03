%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab code to solve a stochastic growth model by PFI
%   Cagliari Summer School, July 2017
%   (c) Antonio Mele
%   (based on G. Hall code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  create the utility function matrices such that, for zero or negative
%  consumption, utility remains a large negative number so that
%  such values will never be chosen as utility maximizing

cons1 = bsxfun(@minus, A_high*kgrid'.^alpha  + delta*kgrid' ,kgrid);
cons2 = bsxfun(@minus, A_low*kgrid'.^alpha  + delta*kgrid' ,kgrid);


cons1(cons1<=0) = NaN;
cons2(cons2<=0) = NaN;

util1 =  log(cons1);
util2 =  log(cons2);

util1(isnan(util1)) = -inf;
util2(isnan(util2)) = -inf;



%%  initialize some variables

v       = zeros(nk,2);   % initial guess for value function:
                            % set to zero for simplicity
tv = zeros(nk,2);   
decis   = zeros(nk,2);   % initial value for policy function
metric  = 10;               % initial value for the convergence metric
iter = 0;
tme = cputime;
[rs,cs] = size(util1);

%%  MAIN LOOP: iterate on Bellman's equation and get the decision rules and
%%  the value function at the optimum


while metric > convcrit;
    
        contv= beta*v*prob'; % continuation value
    
    [tv1,tdecis1]=max(bsxfun(@plus,util1,contv(:,1)) );      
    [tv2,tdecis2]=max(bsxfun(@plus,util2,contv(:,2)) );

    
    tdecis=[tdecis1' tdecis2'];
    
    % Build return vectors
    r1 = zeros(cs,1);
    r2 = zeros(cs,1);
    for i=1:cs
        r1(i) = util1(tdecis1(i),i);
        r2(i) = util2(tdecis2(i),i);
    end
    
    % create matrices Js (see lecture notes)
    g2=sparse(cs,cs);
    g1=sparse(cs,cs);
    for i=1:cs
        g1(i,tdecis1(i))=1;
        g2(i,tdecis2(i))=1;
    end
    % This is the marix P (see lecture notes)
    trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
    
    % Linear algebra step to get the value function associated with P
    tv(:) = ((speye(2*cs) - beta.*trans))\[ r1; r2 ];
    
    jjj
    
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

kgrid = [ (0:ink:maxk)' ];  % capital grid
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state
states   = zeros(n-1,2);
controls = zeros(n-1,2);
[chain,state] = markov(prob,n,s0);
for i = 1:n-1;
    if chain(i) == 1;
        kprime = decis(kmark,1);
        invest = kprime - delta*k;
        cons   = A_high*k^(alpha) - invest;
        kmark = tdecis(kmark,1);
    elseif chain(i) == 2;
        kprime = decis(kmark,2);
        invest = kprime - delta*k;
        cons   = A_low*k^(alpha) - invest;
        kmark = tdecis(kmark,2);
    else
        disp('something is wrong with chain');
    end
    states(i,:) = [ k chain(i) ];
    controls(i,:) = [ cons invest ];
    k = kprime;
end;
