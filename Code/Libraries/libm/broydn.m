% Broydn's method to solve system of nonlinear equations
% INPUT:
%     fname:  name of vector function that should be set to zero
%     xold:   starting value for parameters
%     tolf:   tolerance level
%     iadmat: 1 if Jacobian through admat is available, zero otherwise
%     iprint: 1 if output on iterations is desired
%     additional arguments will be passed to the function
%  OUTPUT:
%     x:      parameter vector that solves the system
%     check:  0 if ok, 1 if there is some problem (no solution found)
%
% The function fname should set the first element of its return vector to
% 1e100 in order to indicate overflow.
function [x, check] = broydn(fname,xold,TOLF,iadmat,iprint,varargin);

  global broydfunc fvec_broydn broydn_in_jacob broydn_step_jacob;

  broydn_in_jacob = 0;
  if(~exist('broydn_step_jacob') | isempty(broydn_step_jacob))
    broydn_step_jacob = 1e-5;
  end;
    
  broydfunc = fname;
  check = 0;

  [nr,nc] = size(xold);
  if(nr>1 & nc>1)
    error('x in broydn must be a vector');
  end;
  x = xold(:);

  MAXITS = 2000; %200;% 50;% 10; %25;% 100; % 5;%
  EPS = 1.0e-7;
  TOLX = EPS;
  STPMX = 100.0;
  TOLMIN = 1.0e-6;
  n=size(x,1);
  
  f=fmin_br(x,varargin{:});
  fx = fvec_broydn;
  if abs(fvec_broydn(1))>=1e100
    error('overflow in function given to broydn at initial vector');
  end;
  test = max(abs(fvec_broydn));
  if (test<TOLF)  %changed from 0.01*TOLF to TOLF;
    return;
  end;
  stpmax=STPMX*max(sqrt(x'*x),n);
  restrt=1;
  since_restrt = 0;
  alam = 1;
  alam2 = 0;
  for its=1:MAXITS
    since_restrt = since_restrt+1;
    txt = sprintf('its:  %d; F''F:  %e\n',its,f);
    if(iprint==1)
      disp(txt(1:end-1));  %without EOL
     
    elseif(iprint==2)
      appendf('broydn_out.txt',txt);
    end;
    % determine whether restart (new computation of Jacobian):
    if(restrt==1 | since_restrt>=3*n) %compute Jacobian at 3n iterations since last restart
      since_restrt = 0;
      alam = 1;  %starting value of backstepping parameter for line search;
      if(iadmat==1)
	xder = deriv(x,eye(n));
	xjac = feval(fname,xder,varargin{:});
	B = getydot(xjac);
      else
	broydn_in_jacob = 1;
	B = jacob(fname,{x,fx},broydn_step_jacob,varargin{:});
	broydn_in_jacob = 0;
      end;
	disp(sprintf('condition number in broydn is %e',cond(B)));
    else
      s = x - xold;
      skip=1;
      w = (fvec_broydn-fvcold) - B*s;
      for i=1:n
	if abs(w(i)) >= EPS*(abs(fvec_broydn(i))+abs(fvcold(i)))
	  skip=0;
	else
	  w(i)=0.0;
	end;
      end;
      if skip==0
	B = B + w*s' / (s'*s);
      end;
    end;
%     if(cond(B)>1e12)
%       error('condition B')
%     end
    p = - B \ fvec_broydn;
    g = B'*fvec_broydn;
    xold = x;
    fvcold = fvec_broydn;
    fold=f;
    [x,f,check,alam2] = lnsrch('fmin_br',xold,fold,g,alam*p,stpmax,varargin{:});
    if(alam2==1)
      alam = 5*alam;
    else
      alam = alam*alam2*2;
    end
    alam = max(min(alam,1),1e-5);
    fx = fvec_broydn;
    test = max(abs(fvec_broydn));
    if (test < TOLF)
      check=0;
      return;
    end;
    if check==1
      if (restrt)
	return;
      else 
	test=0.0;
	den=max(f,0.5*n);
	test = max(abs(g) .* max([abs(x');ones(1,n)])') / den;
	if (test < TOLMIN)
	  return;
	else
	  restrt=1;
	end;
      end;
    else 
      restrt=0;
      test= max(abs(x-xold) ./ max([abs(x');ones(1,n)])');
      if (test < TOLX)
	return;
      end;
    end;
  end;
  check = 1;
  % warning('MXITS exceeded in broydn');
  return;
  
 

