% function [ax,bx,cx] = mnbrak(funcname,ax0,bx0,varargin)
% "brackets" the function 'funcname'
% given two starting value of x (ax and bx), find
% ax < bx < cx such that
% f(bx) < f(ax) and f(bx) < f(cx)
% cf. Numerical Recipes
function [ax,bx,cx,fa,fb,fc] = mnbrak(funcname,ax0,bx0,varargin)
  
  if(iscell(ax0))
    ax = ax0{1};
    fa = ax0{2};
  else
    ax = ax0;
    fa=feval(funcname,ax,varargin{:});
  end
  if(iscell(bx0))
    bx = bx0{1};
    fb = bx0{2};
  else
    bx = bx0;
    fb=feval(funcname,bx,varargin{:});
  end

  if (fb > fa) 
    dum=ax;ax=bx;bx=dum;
    dum=fb;fb=fa;fa=dum;
  end;
  cx=(bx)+1.618034*(bx-ax);
  fc=feval(funcname,cx,varargin{:});
  while (fb > fc) 
    r=(bx-ax)*(fb-fc);
    q=(bx-cx)*(fb-fa);
    u=(bx)-((bx-cx)*q-(bx-ax)*r)/ (2.0*sign2(max(abs(q-r),1.0e-20),q-r));
    ulim=(bx)+100.0*(cx-bx);
    if ((bx-u)*(u-cx) > 0.0) 
      fu=feval(funcname,u,varargin{:});
      if (fu < fc) 
	ax=(bx);
	bx=u;
	fa=(fb);
	fb=fu;
	return;
      elseif (fu > fb) 
	cx=u;
	fc=fu;
	return;
      end;
      u=(cx)+1.618034*(cx-bx);
      fu=feval(funcname,u,varargin{:});
    elseif ((cx-u)*(u-ulim) > 0.0) 
      fu=feval(funcname,u,varargin{:});
      if (fu < fc) 
	bx=cx;cx=u;u=cx+1.618034*(cx-bx);
	fb=fc;fc=fu;fu=feval(funcname,u,varargin{:});
      end;
    elseif ((u-ulim)*(ulim-cx) >= 0.0) 
      u=ulim;
      fu=feval(funcname,u,varargin{:});
    else 
      u=(cx)+1.618034*(cx-bx);
      fu=feval(funcname,u,varargin{:});
    end;
    ax=bx;bx=cx;cx=u;
    fa=fb;fb=fc;fc=fu;
  end;
  




