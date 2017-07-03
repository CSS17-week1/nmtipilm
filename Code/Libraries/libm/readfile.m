function f = readfile(fname)
  if(fileexists(fname))
    f = textread(fname,'%[^\n]');
  else
    f = [];
  end;
