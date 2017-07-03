function f = fmin_br(x,varargin);
  global broydfunc fvec_broydn;
  
  fvec_broydn = feval(broydfunc,x,varargin{:});
  f = 0.5*fvec_broydn'*fvec_broydn;
