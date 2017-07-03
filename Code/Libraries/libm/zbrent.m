function [xmin, fmin] = zbrent(func,x1,x2, tolx,varargin)

  EPS = 1e-8;
  ITMAX = 100;
  if(iscell(x1))
    a=x1{1};
    fa = x1{2}
  else
    a=x1;
    fa = feval(func,a,varargin{:});
  end
  if(iscell(x2))
    b=x2{1};
    fb = x2{2};
  else
    b=x2;
    fb = feval(func,b,varargin{:});
  end
  c=b;

  if ((fa > 0.0 & fb > 0.0) | (fa < 0.0 & fb < 0.0));
    error('Root must be bracketed in zbrent');
    return;
  end;
  fc=fb;
  for (iter=1:ITMAX);
    if ((fb > 0.0 & fc > 0.0) | (fb < 0.0 & fc < 0.0))
      c=a;
      fc=fa;
      d=b-a;
      e=d;
    end;
    if (abs(fc) < abs(fb)) 
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    end;
    tol1=2.0*EPS*abs(b)+0.5*tolx;
    xm=0.5*(c-b);
    if (abs(xm) <= tol1 | fb == 0.0)
      xmin = b;
      fmin = fb;
      return;
    end;
    if (abs(e) >= tol1 & abs(fa) > abs(fb)) 
      s=fb/fa;
      if (a == c) 
	p=2.0*xm*s;
	q=1.0-s;
      else 
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      end;
      if (p > 0.0)
	q = -q;
      end;
      p=abs(p);
      min1=3.0*xm*q-abs(tol1*q);
      min2=abs(e*q);
      if (2.0*p < min(min1,min2)) 
	e=d;
	d=p/q;
      else
	d=xm;
	e=d;
      end
    else 
      d=xm;
      e=d;
    end;
    a=b;
    fa=fb;
    if (abs(d) > tol1)
      b = b + d;
    else
      b = b + sign2(tol1,xm);
    end
    fb = feval(func,b,varargin{:});
  end;
  error('Maximum number of iterations exceeded in zbrent');

  




