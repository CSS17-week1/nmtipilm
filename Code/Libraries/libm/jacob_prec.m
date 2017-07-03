% computes Jacobian finite difference schemes
% Inputs as in jacob
% Outputs:
%   estimate of Jacobian
%   estimate of error
%
% WARNING: the procedure is slow and inefficiently written
% USE ONLY IF COMPUTATION TIME DOES NOT MATTER
function [deriv,err_out] = jacob_prec(func,x_in,h,varargin)

  CON = 1.4;
  CON2 = (CON*CON);
  BIG = 1.0e30;
  NTAB = 10;
  SAFE = 2.0;


  nx = size(x_in,1);
  err_out = [];
  deriv = [];

  if(nargin<3 | h<=0)
    h = 1e-5;
  end;

  for(id=1:nx)
    hh=h;
    x = x_in;
    x(id) = x_in(id)+hh;
    f1 = feval(func,x,varargin{:});

    n = size(f1,1);

    a = nan * ones(n,NTAB,NTAB);
    aux = a;
    
    
    done = zeros(n,1);
  


    x(id) = x_in(id)-hh;
    f0 = feval(func,x,varargin{:});
    f1= f1 - f0;
    f1 = f1 * 1/(2*hh);
    a(:,1,1)= f1;
    err = BIG * ones(n,1);
    errt = err;
    ans = err;
    for (i=2:NTAB) 
      hh = hh / CON;
      
      x(id) = x_in(id)+hh;
      f1 = feval(func,x,varargin{:});
      x(id) = x_in(id)-hh;
      f0 = feval(func,x,varargin{:});
      f1= f1 - f0;
      f1 = f1 * 1/(2*hh);
      a(:,1,i)= f1;
      
      fac=CON2;
      for (j=2:i)
	aux = a(:,j-1,i);
	aux = aux * fac;
	aux = aux - a(:,j-1,i-1);
	aux = aux/(fac-1.0);
	a(:,j,i) = aux;
	fac=CON2*fac;
	for(l=1:n)
	  if(done(l)==0)
	    errt(l)=max(abs(a(l,j,i)-a(l,j-1,i)), ...
			 abs(a(l,j,i)-a(l,j-1,i-1)));
	    if (errt(l) <= err(l)) 
	      err(l)=errt(l);
	      ans(l)=a(l,j,i);
	    end;
	  end;
	end;
      end;
      for(l=1:n)
	if (abs(a(l,i,i)-a(l,i-1,i-1)) >= SAFE*(err(l)))
	  done(l) = 1;
	end;
      end;
    end;

    err_out = [err_out err];
    deriv = [deriv ans];
  end;
