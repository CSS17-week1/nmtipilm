function B = funbas2(fspace,x,order);
  if(iscell(fspace.bastype))
    if(nargin==2)
      B = funbas(fspace,x);
    else
      B = funbas(fspace,x,order);
    end;
  else
    switch fspace.bastype
      case 'complpoly'
	itensor = 0;
      case 'tenspoly'
	itensor = 1;
      otherwise
	error('wrong bastype in funeval2');
    end;
    if(nargin>2 & any(order~=0))
      error('order must be zero for bastypes complpoly or tenspoly');
    end;
    c = zeros(nparams(fspace.n-1,itensor));
%         [unused,B] = chebeval(x,c,fspace.a,fspace.b,fspace.n-1);

    [unused,B] = fastchebeval(x,c,fspace.a,fspace.b,fspace.n-1,itensor);
  end;












