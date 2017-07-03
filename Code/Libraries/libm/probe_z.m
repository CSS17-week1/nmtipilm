function z = probe_z(x)
  global probe_z_ncalls;
  probe_z_ncalls = probe_z_ncalls +1;

%  z = x.^4 - 16;
%  return;

%  y = min(max(x,-100),100);
%  z = exp(y)-exp(9);
%  return;

  if(x<-0.001)
    z = -1;
  elseif(x>0.001)
    z = 1;
  else
    z = (x*1000)^3;
  end;
  z = -z;