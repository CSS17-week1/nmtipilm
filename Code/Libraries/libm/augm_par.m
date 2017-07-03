function par2 = augm_par(par,apprord,apprordnew,itensor);
  
  n1 = apprord(1)-1;
  n2 = apprord(2)-1;
  n1new = apprordnew(1)-1;
  n2new = apprordnew(2)-1;
  npar = nparams(apprord-1,itensor);
  if(npar~=size(par,1) | size(par,2)~=1)
    error('input array of wrong size');
  end;
  npar2 = nparams(apprordnew-1,itensor);
  if(itensor==0)
    if(n1~=n2 | n1new~=n2new);
      error('approrder must be same in both dimensions for compl.pol');
    end;
%    par2 = [par;zeros(npar2-npar,1)];
    par2 = zeros(npar2,1);
    for i=0:n1;
      for j=0:n1-i;
	par2(polyindx([i,j],[n1new,n2new],itensor)) ...
	    = par(polyindx([i,j],[n1,n2],itensor));
      end;
    end;
  else
    par2 = zeros(npar2,1);
    for i=0:n1;
      for j=0:n2;
	par2(polyindx([i,j],[n1new,n2new],itensor)) ...
	    = par(polyindx([i,j],[n1,n2],itensor));
      end;
    end;
  end;
