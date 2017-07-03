% unconstrained optimization using BFGS with ldl factorization
% INPUT:
%     func:   name of vector function that should be set to zero
%     xold:   starting value for parameters
%     gtol:   tolerance level for gradient
%     ideriv: how to calculate the gradient
%             0 for forward differences (jacob.m)
%             1 for ADMAT
%             2 for higher order difference scheme (jacob_prec.m)
%     iprint: 1 if output on iterations is desired
%     additional arguments will be passed to the function
%  OUTPUT:
%     p:      optimal parameter vector found
%     fret:   best function value found 
%     g:      gradient at p
%     iter:   number of iterations needed
%
% If some range error (overflow) occurs in a function evaluation,
% return HUGE number (e.g. 1e100)
% Michael Reiter, last change 5.1.2004
function [p,fret,g,iter] = bfgs(func,xold, gtol, ideriv, iprint, varargin)

  ITMAX = 200;
  EPS = 3.0e-8;
  TOLX = (4*EPS);
  STPMX = 100.0;
  
  p = xold;
  n = size(p,1);
  dg = zeros(n,1);
  g = zeros(n,1);
  hdg = zeros(n,1);
  hess_ldl = eye(n);
  pnew = zeros(n,1);
  xi = zeros(n,1);
  fp= feval(func,p,varargin{:});
  if(ideriv==0)
    g = jacob(func,p,1e-5,varargin{:});
  elseif(ideriv==1)
    pder = deriv(p,eye(size(p,1)));
    pjac = feval(func,pder,varargin{:});
    g = getydot(pjac);
  elseif(ideriv==2)
    g = jacob_prec(func,p,1e-5,varargin{:});
  else
    error('inadmissible value for ideriv');
  end;

    
  g = g(:);  %make column vector;
  xi = -g; %steepest descent
  sum = p'*p;
  stpmax=STPMX*max(sqrt(sum),n);
  for (its=1:ITMAX) 
    iter=its;
    [pnew,fret,check] = lnsrch(func,p,fp,g,xi,stpmax,varargin{:});
    fp = fret;
    xi=pnew-p; %(p is x!!) change in x
    p=pnew;
    
    test = max(abs(xi)./max(abs(p),1));
    if (test < TOLX) 
      return;
    end
    dg=g;
    if(ideriv==0)
      g = jacob(func,p,1e-5,varargin{:});
    elseif(ideriv==1)
      pder = deriv(p,eye(size(p,1)));
      pjac = feval(func,pder,varargin{:});
      g = getydot(pjac);
    elseif(ideriv==2)
      g = jacob_prec(func,p,1e-5,varargin{:});
    else  
      error('inadmissible value for ideriv');
    end;
    g = g(:);  %make column vector;
    den=max(abs(fret),1.0);
    test=max(abs(g).*max(abs(p),1.0)/den);%elasticity df/dx
    if iprint
        disp(sprintf('iter: %d; func: %e; max grad:  %e; check=%d',its,fret,test,check));
    end;
    if (test < gtol) 
      return;
    end
    dg=g-dg;

    fac =   dg'*xi;
    sumdg = dg'*dg;
    sumxi = xi'*xi;
    if (fac*fac > EPS*sumdg*sumxi) 
      % update; the ldl factorization used guarantees
      % a pos.def. hessian in any case;
      hessxi = ldlmul(hess_ldl,xi);
      fah = xi'*hessxi;
      hess_ldl = ldlupd(hess_ldl,dg,1/fac);
      hess_ldl = ldlupd(hess_ldl,hessxi,-1./fah);
    end

    xi = - ldlsolve(hess_ldl,g);
  end
  warning('too many iterations in bfgs');

