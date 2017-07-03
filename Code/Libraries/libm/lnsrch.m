function [x,f,check,alam] = lnsrch(fname,xold,fold,g,p,stpmax,varargin);
  
  n = size(xold,1);
  ALF = 1.0e-4;
  TOLX = 1.0e-7;
  check=0;
  sum = p'*p;
  sum=sqrt(sum);
  if (sum > stpmax)
    p = p * stpmax/sum;
  end;
  slope = g'*p;
  test=max(abs(p) ./ max([abs(xold');ones(1,n)])');
  % min() introduced by MR 5.4.2006:
  % otherwise, alamin can be >1 !!
  alamin=min(0.1,TOLX/test);  
  alam=1.0;
  for its=1:100000
    x=xold + alam*p;
    try
      f=feval(fname,x,varargin{:});
    catch
      f = 1e100;  %something outrageous;
    end;
    if (alam < alamin) 
      x = xold;
      check=1;
      return;
    else
      if f <= fold+ALF*alam*slope
	return;
      else 
	if alam == 1.0
	  tmplam = -slope/(2.0*(f-fold-slope));
	else 
	  rhs1 = f-fold-alam*slope;
	  rhs2=f2-fold2-alam2*slope;
	  a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	  b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	  if (a == 0.0)
	    tmplam = -slope/(2.0*b);
	  else 
	    disc=b*b-3.0*a*slope;
	    if (disc<0.0)
	      error('Roundoff problem in lnsrch');
	    else
	      tmplam=(-b+sqrt(disc))/(3.0*a);
	    end;
	  end;
	  if (tmplam>0.5*alam)
	    tmplam=0.5*alam;
	  end;
	end;
      end;
    end;
    alam2=alam;
    f2 = f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
  end;
  


