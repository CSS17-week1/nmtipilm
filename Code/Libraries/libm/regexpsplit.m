function x = regexpsplit(str,patt)
  [i0,i1] = regexp(str,patt);
  n = length(i0);
  x = cell(n+1,1);
  len = length(str);
  if(n==0)
    x{1} = str;
  else
    x{1} = str(1:i0(1)-1);
    for i=2:n
      x{i} = str(i1(i-1)+1:i0(i)-1);
    end
    x{n+1} = str(i1(n)+1:len);
  end
