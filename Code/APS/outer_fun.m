function c = outer_fun(w,h,delta,pi)

% w = future payoffs
% delta = discount factor
% H,c ineqaulities
% h = current search subgradient
% sgp = current payoff


c_tmp = h(1) * ((1-delta)*pi(1)+delta*w(1)) + h(2) * ((1-delta)*pi(2)+delta*w(2));

c = - c_tmp;