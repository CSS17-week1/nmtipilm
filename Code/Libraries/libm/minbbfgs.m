% Minimize a function f(x) s.t. simple bound constraints x>=xlow
% Theory: Bertsekas, "Nonlinear Programming", Section 1.5
% Inputs:
%   func:        objective function 
%   x0:          initial guess for x
%   xlow:        lower bound for x
%   TOLG:        accept point if abs(gradient)<TOLG
%   iAdmat:      if true, compute gradient by admat, otherwise by forw.diff.
%   imonitor:    if true, print info at each iteration 
%   varargin:    any additional parameters are passed on to f
% Outputs:
%   x:           minimum in case of sucess, otherwise best estimate
%   info:        0 (if success)
%                1 (if it fails during line search)
%                2 (if maximum iteration number is exceeded)
%   slack:       complementary slackness term:
%                  max(abs((x-xlow).*g))
%  Michael Reiter, October 2004
%  Last update:  October 2005
function [x,info,slack] = minbbfgs(func,x0, xlow, TOLG, iAdmat, imonitor, varargin)
  global B;  % make global to check estimate after convergence;
  % we say constraint binds if (x(i)-xlow(i))<eps:
  eps = 1e-8; %typical x-values should be around 1!;
  eps_update = 3.0e-8;  %cf. Numerical recipes, bfgs algorithm
  sigma = 1e-4; %parameter in linear search algorithms
  maxits = 2000; %max. number of iterations

  if(size(x0,2)~=1 | size(xlow,2)~=1)
    error('starting guess and lower bounds in minbounds must be column vector');
  end;
  n = size(x0,1);
  if(n~=size(xlow,1))
    error('x0 and xlow must have same size in minbounds');
  end;
  x = max(x0,xlow);

  f = feval(func,x,varargin{:});
  if(f>=1e100)
    error('initial function value in minbbfgs>=1e100; inadmissible input');
  end;

  info = 0; %we assume for the moment that the result is ok (until bad things happen);
  B = eye(n); %initialize estimate of H;
  if(iAdmat)
    xder = deriv(x,eye(n));
    tmp = feval(func,xder,varargin{:});
    g = getydot(tmp)';
  else
    g = jacob(func,x,1e-5);
    g = g(:);  % make it column vector
  end;
  for its=1:maxits
    slack = max(abs((x-xlow).*g)); 
    if(all(abs(g)<TOLG | (abs(x)<=xlow+eps & g>0)))%close enough to optimum;
      return;
    end;

    % index of elements with binding constraint:
    I = x<xlow+eps & g>0;
    nconstr = sum(I); %number of constraints binding;
    NotI = ~I;  %index of elements with nonbinding constraint;
    if(imonitor)
      disp(sprintf('Iter %d: f=%e; #Constr=%d; slack=%e',...
		   its,f,nconstr,slack));
    end;

    if(nconstr==n) %all variables constrained
      return;
    end;
    if(nconstr==0)
      p = -B\g;
    else
      p = zeros(n,1);
      p(NotI) = -B(NotI,NotI)\g(NotI); %Newton step for the free parameters
      p(I) = - g(I)./diag(B(I,I));  %points into infeasible set;
    end;


    if( max(abs(p(NotI)))<1e-8)%Newton step tiny; should be close to optimum;
      return;
    end;

    xold = x;
    alpha = 1;
    arm_ok=0;
    while(alpha>1e-9)
      dx = alpha*p;
      xtry = max(x+dx,xlow);
      ftry = feval(func,xtry,varargin{:});
      df_expect = sum(g.*(I.*(xtry-x) + NotI.*dx));
      if(ftry-f < sigma*df_expect) %Armijo rule satisfied
    	x = xtry;
    	f = ftry;
    	arm_ok=1;
	break;
      else %otherwise step back
	alpha = 0.5*alpha;
      end;
    end;
    if(arm_ok==0) %no acceptable step size found;
      info = 1;
      return;
    end;

    gold = g;

    if(iAdmat)
      xder = deriv(x,eye(n));
      tmp = feval(func,xder,varargin{:});
      g = getydot(tmp)';
    else
      g = jacob(func,x,1e-5);
      g = g(:);
    end;

    % update of Hessian estimate:
    xi = x-xold;
    dg = g - gold;
    fac =   dg'*xi;
    sumdg = dg'*dg;
    sumxi = xi'*xi;
    if (fac*fac > eps_update*sumdg*sumxi) 
      % update; the ldl factorization used guarantees
      % a pos.def. hessian in any case;
      Bxi = B*xi;
      C = dg*dg'/(dg'*xi) - Bxi*Bxi'/(xi'*Bxi);
      B = B + C;
    end
  end;



  info = 2; %maxits exceeded;
