% Minimize a function f(x) s.t. simple bound constraints x>=xlow
% Theory: Bertsekas, "Nonlinear Programming", Section 1.5
% Inputs:
%   x0:          initial guess for x
%   xlow:        lower bound for x
%   func:        function that gives
%                   f, if called with 1 output argument
%                   [gradient, hessian], if called with 2 output arguments
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
%  Last update:  -
function [x,info,slack] = minbounds(func, x0, xlow, TOLG, imonitor, varargin)

  eps = 1e-7; %typical x-values should be around 1!;
  sigma = 1e-4;
  maxits = 2000;

  if(size(x0,2)~=1 | size(xlow,2)~=1)
    error('starting guess and lower bounds in minbounds must be column vector');
  end;
  n = size(x0,1);
  if(n~=size(xlow,1))
    error('x0 and xlow must have same size in minbounds');
  end;
  x = max(x0,xlow);

  f = feval(func,x,varargin{:});

  info = 0; %we assume for the moment that the result is ok (until bad things happen);
  for its=1:maxits
    %2 output arguments to func mean that we need gradients:
    [f,g,H] = feval(func,x,varargin{:});  
    g = g(:);
    if(all(abs(g)<TOLG | (abs(x)<=xlow+eps & g>0)))%close enough to optimum;
      return;
    end;


    slack = max(abs((x-xlow).*g)); 
    I = x<xlow+eps & g>0;
    nconstr = sum(I);
    NotI = ~I;
    if(imonitor)
      disp(sprintf('Iter %d: f=%e; #Constr=%d; slack=%e',...
		   its,f,nconstr,slack));
    end;

    if(nconstr==n) %all variables constrained
      return;
    end;
    if(nconstr==0)
      p = -H\g;
    else
      p = zeros(n,1); % enforce column vector;
      p(NotI) = -H(NotI,NotI)\g(NotI);
      p(I) = - g(I)./diag(H(I,I));
    end;


    if( max(abs(p(NotI)))<1e-8)%Newton step tiny; should be close to optimum;
      return;
    end;

    alpha = 1;
    arm_ok=0;
    while(alpha>1e-9)
      dx = alpha*p;
      xtry = max(x+dx,xlow);
      ftry = feval(func,xtry,varargin{:});
      df_expect = sigma*sum(g.*(I.*(xtry-x) + NotI.*dx));
      if(ftry-f < df_expect)
    	x = xtry;
    	f = ftry;
    	arm_ok=1;
	break;
      end;
      alpha = 0.5*alpha;
    end;
    if(arm_ok==0)%not acceptable step size found;
      info = 1;
      return;
    end;
  end;

  info = 2; %maxits exceeded;
