function ex = fileexists(fn)
  ex = 0;
  fid = fopen(fn);
  if(fid>=0)
    ex = 1;
    fclose(fid);
  end;
