function ds = dir_separator
  m = matlabroot;
  for i=1:length(m)
    if(m(i) == '/')
      ds = '/';
      return;
    elseif(m(i) == '\')
      ds = '\';
      return;
    end
  end
  error('cannot determine directory separator');
