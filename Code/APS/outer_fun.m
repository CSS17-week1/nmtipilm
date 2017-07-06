<<<<<<< HEAD
<<<<<<< HEAD
function c = outer_fun(w,h,delta,pi)

% w = future payoffs
% delta = discount factor
% H,c ineqaulities
% h = current search subgradient
% sgp = current payoff


c_tmp = h(1) * ((1-delta)*pi(1)+delta*w(1)) + h(2) * ((1-delta)*pi(2)+delta*w(2));

=======
function c = outer_fun(w,h,delta,pi)

% w = future payoffs
% delta = discount factor
% H,c ineqaulities
% h = current search subgradient
% sgp = current payoff


c_tmp = h(1) * ((1-delta)*pi(1)+delta*w(1)) + h(2) * ((1-delta)*pi(2)+delta*w(2));

>>>>>>> 0b5165c8e0cbd0e4f1dde7640d267055799a76ba
=======
function c = outer_fun(w,h,delta,pi)

% w = future payoffs
% delta = discount factor
% H,c ineqaulities
% h = current search subgradient
% sgp = current payoff


c_tmp = h(1) * ((1-delta)*pi(1)+delta*w(1)) + h(2) * ((1-delta)*pi(2)+delta*w(2));

>>>>>>> 0b5165c8e0cbd0e4f1dde7640d267055799a76ba
c = - c_tmp;