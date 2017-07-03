function [y,Ch] = funevalpoly(c,fspace,x);
  itensor = 0;
  if(nargout==2)
    [y,Ch] = fastchebeval(x,c,fspace.a,fspace.b,fspace.n-1,itensor);
  else
    y = fastchebeval(x,c,fspace.a,fspace.b,fspace.n-1,itensor);
  end;
end;












