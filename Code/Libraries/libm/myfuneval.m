function [y,Ch] = funeval2(c,fspace,x);
  if(iscell(fspace.bastype))
    y = funeval(c,fspace,x);
  else
    switch fspace.bastype
      case 'complpoly'
	itensor = 0;
      case 'tenspoly'
	itensor = 1;
      otherwise
	error('wrong bastype in funeval2');
    end;
    if(nargout==2)
      [y,Ch] = fastchebeval(x,c,fspace.a,fspace.b,fspace.n-1,itensor);
    else
      y = fastchebeval(x,c,fspace.a,fspace.b,fspace.n-1,itensor);
    end;
  end;












