% function [fmin, xmin] = brent(func,ax, bx, cx, tol);
% finds minimum of function 'func',
% arguments ax,bx,and cx are the output of mnbrak
% tol is the precision of x on return
% cf. Numerical Recipes
function [xmin, fmin] = brent(func,ax, bx, cx, tol,varargin);

  ITMAX = 100;
  CGOLD = 0.3819660;
  ZEPS = 1.0e-10;

  e=0.0;

  a=min(ax,cx);
  b=max(ax,cx);
  if(iscell(bx))
    fx = bx{2};
    bx = bx{1};
  else
    fx = feval(func,bx,varargin{:});
  end
  v=bx;
  w=v;
  x=w;
  fv=fx;
  fw=fv;
  for iter=1:ITMAX
    xm=0.5*(a+b);
    tol1=tol*abs(x)+ZEPS;
    tol2=2.0*(tol1);
    if (abs(x-xm) <= (tol2-0.5*(b-a))) 
      xmin=x;
      fmin = fx;
      return;
    end;
    if (abs(e) > tol1) 
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0)
	p = -p;
      end;
      q=abs(q);
      etemp=e;
      e=d;
      if (abs(p) >= abs(0.5*q*etemp) | p <= q*(a-x) | p >= q*(b-x))
	if x >= xm
	  e=a-x;
	else
	  e= b-x;
	end;
	d=CGOLD*e;
      else 
	d=p/q;
	u=x+d;
	if (u-a < tol2 | b-u < tol2)
	  d=sign2(tol1,xm-x);
	end;
      end;
    else 
      if x >= xm
	e=a-x;
      else
	e= b-x;
      end;
      d=CGOLD*e;
    end;
    if abs(d) >= tol1
      u= x+d;
    else
      u= x+sign2(tol1,d);
    end;
    fu = feval(func,u,varargin{:});
    if (fu <= fx) 
      if (u >= x)
	a=x;
      else
	b=x;
      end;
      v=w;
      w=x;
      x=u;
      fv=fw;
      fw=fx;
      fx=fu;
    else 
      if (u < x)
	a=u;
      else
	b=u;
      end;
      if (fu <= fw | w == x) 
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      elseif (fu <= fv | v == x | v == w) 
	v=u;
	fv=fu;
      end;
    end;
  end;
  error('Too many iterations in brent');
  xmin=x;
  fmin = fx;
  return;


