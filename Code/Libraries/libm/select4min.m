function f = select4min(lambda,func,x,delta,varargin)
  f = feval(func,x+lambda*delta,varargin{:});
