% PROVIDE COMMENTS ON THOSE PLACES WHERE YOU FIND '??'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Quadratic-Hill-Climbing-Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:   ff:  function to minimize
%          p:   starting value
%Output:   optimal x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x0,check] = qhc(ff,x,tolg,iprint,varargin)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %          Parameters:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  R = 0.02;
  maxit = 200;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %          Initialize some values:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n = size(x,1);
  x0 = x;  
  x1 = x;
  fw1 = ff(x0,varargin{:});
  check = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %           Iterations:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iter=1:maxit
    x0 = x1;
    fw0 = fw1;
    [fw0,g,H] = ff(x0,varargin{:});
    g = g(:);
    relgrad = (g .* max(x0,1)) / max(fw0,1);    % gradient, as elasticitiy
    lambda = min(eig(H));
    if all(abs(relgrad)<tolg)
      if lambda > 0
	return;
      else
	% COMMENT: WHAT KIND OF POINT HAVE WE REACHED HERE??
	[evec,eval] = eig(H);
	eval = diag(eval);
	[lambda,il] = min(eval);
	d = evec(:,il);
	% COMMENT: WHY DO WE CHOOSE THE SEARCH DIRECTION D??
	d = d * max(1,norm(x0,2)) / norm(d,2);
      end
    else
      if(lambda>1e-8)
	a = 0;
      else
	a = abs(lambda) + R*norm(g,2);
      end;
      S = H + (a.*eye(n));      
      % Quasi-Newton search direction:
      d = - S\g;
    end

    step = 1;
    search_ok = 0;
    while(step>1e-8)
      x1 = x0 + step*d;
      fw1 = ff(x1,varargin{:});
      if(step==1)
	fwq = fw0 + d'*g + 0.5*d'*H*d;
	% COMMENT: WHAT IS THE MEANING OF Z??
	z = (fw1 - fw0)/(fwq - fw0);
	% COMMENT: WHY DO WE AJUST R IN THE FOLLOWING WAY??
	if z <= 0;             
	  R = 4*R;
	elseif (0<z) & (z<=0.8);
	  R = (4 - 3.6*z) * R;    
	elseif (0.8<z) & (z<=1.2);
	  R = 0.4*R;              
	elseif (1.2<z) & (z<=2); 
	  R = (-3.2 + 3.6*z) * R; 
	else %z > 2;  
          R = 4*R; 
        end;
      end
      % COMMENT: WHAT IS THE MEANING OF THE TERM step*(d'*g)??
      % (the following criterion, rather than fw1<fw0,
	% avoids steps with smaller and smaller improvements)
      if (fw1 < fw0 - 1e-4*step*(d'*g))
	search_ok = 1;
	break;
      end;
      % stepping back:
      step = 0.25*step;
    end
    if(iprint)
      disp(sprintf('iter=%d,lambda=%g, R = %g, f = %g, g = %e',...
		   iter,lambda, R, fw0, max(abs(relgrad))));
      % x1'
    end;
    if(search_ok==0)
      check = 1;
      return;
    end;
  end;
  check = 2;
  return;

