function G = funnode2(fspace);
  if(iscell(fspace.bastype))
    G = funnode(fspace);
  elseif ( strcmp(fspace.bastype,'complpoly') | strcmp(fspace.bastype,'tenspoly'))
    fspace2 = fundefn('cheb',fspace.n,fspace.a,fspace.b);
    G = funnode(fspace2);
  else
    error('wrong bastype in funeval2');
  end;












