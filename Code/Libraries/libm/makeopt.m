function opt = solveopt(fname)
  switch(fname)
    case 'solve1d'
      opt = struct('ftol',[],'xtol',[],'bounds',[]);
    otherwise
      error('unknown optimization routine');
  end;
