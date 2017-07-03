function [p,fret,iter] =  powell(func, p, xi,ftol,maxstep,varargin)
  
  n = length(p);
  if(isempty(xi))
    xi = eye(n);
  end

  maxiter = 200;
  fret=feval(func,p,varargin{:});
  pt = p;
  for iter=1:maxiter
    fp=fret;
    ibig=0;
    del=0.0;
    for (i=1:n) 
      xit=xi(:,i);
      fptt=fret;
      step = maxstep/max(abs(xit));
      %[lambda,fret] = minbnd(@select4min,-step,step,1e-6,func,p,xit,varargin{:});
      %[lambda,fret] = fminbnd(@select4min,-step,step,[],func,p,xit,varargin{:});
      [ax,bx,cx,fa,fb] = braksimple(@select4min,-step,step,func,p,xit,varargin{:});
      tol = 1e-5;
      [lambda,fret] = brent(@select4min,ax, {bx,fb}, cx, tol,func,p,xit,varargin{:});

      p = p+lambda*xit;
      % disp(sprintf('fret = %e',fret));
      if (abs(fptt-fret) > del) 
	del=abs(fptt-fret);
	ibig=i;
      end
    end
    if (2.0*abs(fp-fret) <= ftol*(abs(fp)+abs(fret))) 
    % if (2.0*abs(fp-fret) <= ftol)
      return;
    end
    if (iter == maxiter) 
      error('powell exceeding maximum iterations');
    end
    ptt=2.0*p-pt;
    xit=p-pt;
    pt=p;
    fptt=feval(func,ptt,varargin{:});
    if (fptt < fp) 
      t=2.0*(fp-2.0*fret+fptt)*(fp-fret-del)^2-del*(fp-fptt)^2;
      if (t < 0.0) 
	step = maxstep/max(abs(xit));
	%[lambda,fret] = minbnd(@select4min,-step,step,1e-6,func,p,xit,varargin{:});
	[lambda,fret] = fminbnd(@select4min,-step,step,[],func,p,xit,varargin{:});
	p = p+lambda*xit;
	xi(:,ibig)=xi(:,n);
	xi(:,n)=xit;
      end
    end
  end

