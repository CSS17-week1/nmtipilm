function z = select4zero(xi,func,x,ix,iy,varargin)
  x(ix) = xi;
  zz = feval(func,x,varargin{:});
  z = zz(iy);