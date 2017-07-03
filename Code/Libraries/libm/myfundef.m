function fspace = myfundef(bastype,n,a,b)
  if(strcmp(bastype,'complpoly') | strcmp(bastype,'tenspoly'))
    d = length(a);
    fspace = struct('d',d,'n',n,'a',a,'b',b,'bastype',bastype);
  else
    fspace = fundefn(bastype,n,a,b);
  end;
