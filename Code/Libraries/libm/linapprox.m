% nPerDim: number of grid points per dim
% a: lower bounds
% b: upper bounds
% ilog, logshift:
%     if ilog(i)=1, grid is uniform in log(x(i)+logshift(i)),
%     otherwise uniform in x
% iMultiLin: if 1, multilinear interpolation, otherwise
%      simplicial linear interpolation
function fspace = linapprox(nPerDim,a,b,ilog,logshift,iMultiLin)
  d = length(nPerDim);
  nmax = max(nPerDim);
  grid = zeros(d,nmax);
  gridinvdist = zeros(d,nmax-1);


  nPerDimCumul = zeros(d,1);
  nPerDimCumul(d) = 1;
  for i=d-1:-1:1
    nPerDimCumul(i) = nPerDim(i+1) * nPerDimCumul(i+1);
  end;

  for i=1:d
    if(ilog(i)==1)
      x = exp(linspace(log(a(i)+logshift(i)),log(b(i)+logshift(i)), ...
		       nPerDim(i))) - logshift(i);
      x(1) = a(i);
      x(nPerDim(i)) = b(i);
    else
      x = linspace(a(i),b(i),nPerDim(i));
    end;
    grid(i,1:nPerDim(i)) = x;
    gridinvdist(i,1:nPerDim(i)-1) = 1./diff(x);
  end;

  if(iMultiLin==1)
    iML = 1;
  else
    iML = 0;
  end;
  fspace = struct('d',d,'n',nPerDim,'nCumul',nPerDimCumul,...
		  'a',a,'b',b,'grid',grid,...
		  'gridInvDist',gridinvdist,'iML',iML);

