n = 200;

Q = eye(n);


midp = ((1:n)-0.5)/n;

H = [ones(1,n);midp];
m = [1;0.8];

%
% Augm. Lagrangian is 
% 0.5*f'Qf + lambda'H*f + 0.5*c*(H*f-m)'(H*f-m)
%
%


lambda = zeros(size(H,1),1);
c = 10000;


for i=1:10
  A = Q + c*H'*H;
  b = H'*lambda - c*H'*m;
  [f,val,g,its] = qpsimple(A,b,ones(n,1));
  ConstrViol = H*f-m;
  disp(its);
  disp(ConstrViol');
  if(all(abs(ConstrViol)<1e-12))
    break;
  end;
  lambda = lambda + c*ConstrViol;
  disp(lambda');
end;
