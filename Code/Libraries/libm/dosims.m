function [G1,impact,eu,abseig] = dosims(funcequ,stst,ix,ixlag,ieps,ieta)

  ntotal = max([max(ix),max(ixlag),max(ieps),max(ieta)]);

  stst2 = zeros(ntotal,1);
  stst2([ixlag ix]) = repmat(stst,2,1);
  resid = feval(funcequ,stst2);
  check = max(abs(resid));
  if(check>1e-10)
    disp('Maximum residual at stst:');
    disp(check);
    pause;   % ALWAYS CHECK THAT THIS IS VERY SMALL!!;
  end;

  jac = jacob(funcequ,stst2,1e-6);

  % SAME AS IN ALL OTHER APPLICATIONS:
  g0 = -jac(:,ix);
  g1 = jac(:,ixlag);
  c = zeros(size(g0,1),1);
  Psi = jac(:,ieps);
  Pi = jac(:,ieta);
  div = 1+1e-6

  [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,Psi,Pi,div);
  if(~isempty(gev))
    abseig = sort(abs(gev(:,2)./gev(:,1)));
  else
    abseig = [];
  end;



