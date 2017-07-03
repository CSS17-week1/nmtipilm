function [xmin,f,g,its] = qpsimple(Q,b,xstart)
  
  eps = 1e-4;
  sigma = 1e-4;

  n = size(Q,1);

  x = xstart;
  f = 0.5*x'*Q*x + b'*x;

  isuccess = 0;

  for its=1:1000
    g = Q*x + b;
    I = x<eps & g>0;
    II = min(I*ones(1,n) + ones(n,1)*I',1) -diag(I,0);
    H = Q-II.*Q;
    p = H\g;

    KT = abs(g)<1e-7 | (abs(x)<1e-7 & g>0);
    if(all(KT))
      isuccess = 1;
      break;
    end;

    alpha = 1;

    arm_ok=0;
    for j=1:40
      xtry = x-alpha*p;
      xtry = max(xtry,0);
      ftry = 0.5*xtry'*Q*xtry + b'*xtry;
      crit = sigma*g'*(alpha*p.*(1-I) + I.*(x-xtry));
      if f-ftry > crit;
	x = xtry;
	f = ftry;
	arm_ok=1;
	break;
      end;
      alpha = 0.5*alpha;
    end;
    if(arm_ok==0)
      warning('armijo rule never accepted');
    end;
%    disp(f);
  end;

  xmin = x;
  f = 0.5*x'*Q*x + b'*x;
  g = Q*x + b;

  if(isuccess==0)
    disp(sprintf('Warning: x*g=%e',x'*g));
  end;

