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
    u=(cx)+1.618034*(cx-bx);
    fu=feval(funcname,u,varargin{:});
    ax=bx;bx=cx;cx=u;
    fa=fb;fb=fc;fc=fu;
  end;
  




