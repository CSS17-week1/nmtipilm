function [p,fret,g,iter] = dfpmin(p, gtol, func, varargin)

  ITMAX = 200;
  EPS = 3.0e-8;
  TOLX = (4*EPS);
  STPMX = 100.0;
  
  n = size(p,1);
  dg = zeros(n,1);
  g = zeros(n,1);
  hdg = zeros(n,1);
  hessin = eye(n);
  pnew = zeros(n,1);
  xi = zeros(n,1);
  fp= feval(func,p,varargin{:});
  g = jacob(func,p,1e-5,varargin{:});
  g = g';
  xi = -g;
  sum = p'*p;
  stpmax=STPMX*max(sqrt(sum),n);
  for (its=1:ITMAX) 
    iter=its;
    [pnew,fret,check] = lnsrch(func,p,fp,g,xi,stpmax,varargin{:});

    fp = fret;
    xi=pnew-p;
    p=pnew;
    
    test=0.0;
    for (i=1:n) 
      temp=abs(xi(i))/max(abs(p(i)),1.0);
      if (temp > test) 
	test=temp;
      end
    end
    if (test < TOLX) 
      return;
    end
    dg=g;
    g = jacob(func,p,1e-5,varargin{:});
    g = g';
    test=0.0;
    den=max(fret,1.0);
    for (i=1:n) 
      temp=abs(g(i))*max(abs(p(i)),1.0)/den;
      if (temp > test) 
	test=temp;
      end
    end
    if (test < gtol) 
      return;
    end
    dg=g-dg;
    
    hdg = hessin*dg;
    
    fac = dg'*xi;
    fae = dg'*hdg;
    sumdg = dg'*dg;
    sumxi = xi'*xi;
    if (fac*fac > EPS*sumdg*sumxi) 
      fac=1.0/fac;
      fad=1.0/fae;
      dg=fac*xi-fad*hdg;
      for (i=1:n) 
	for (j=1:n) 
	  hessin(i,j) = hessin(i,j) + fac*xi(i)*xi(j) ...
	  -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
	end
      end
    end
    xi = -hessin*g;
  end
  warning('too many iterations in dfpmin');

