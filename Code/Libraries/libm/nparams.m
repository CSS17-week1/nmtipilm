function n = nparams(polorder,itensor)
  if(itensor==1)
    n = prod(polorder+1);
  else
    if(size(polorder,1)==1)
      polorder = polorder';
    end;
    dim = size(polorder,1);
    n = prod(polorder + (1:dim)') / factorial(dim);
  end;
