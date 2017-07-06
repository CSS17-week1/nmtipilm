<<<<<<< HEAD
function [w,c,flag]=inner_fun(h,nv,v,pi,w_min,pi_star,delta)

% Choosing future values for each action
% v are the vertices (reduced before)
% nv number of vertices
% h search subgradient
% current payoff

alpha = 0:.01:1;
na = length(alpha);
T = (nv-1)*na;

c_tmp = zeros(T,1);
w_tmp = zeros(T,2);
flag_tmp = zeros(T,1);

j=1;
for aa=1:na
    for ih=1:nv-1
        w1 = (1-alpha(aa))*v(ih,1)+alpha(aa)*v(ih+1,1);
        w2 = (1-alpha(aa))*v(ih,2)+alpha(aa)*v(ih+1,2);
        
        ic_1 = (1-delta)*pi(1)+delta*w1-(1-delta)*pi_star(1)-delta*w_min(1);
        ic_2 = (1-delta)*pi(2)+delta*w2-(1-delta)*pi_star(2)-delta*w_min(2);
        if ic_1>=0 & ic_2>=0
            c_tmp(j)    = h(1) * ((1-delta)*pi(1)+delta*w1) + h(2) * ((1-delta)*pi(2)+delta*w2);
            w_tmp(j,1)  = (1-delta)*pi(1)+delta*w1;
            w_tmp(j,2)  = (1-delta)*pi(2)+delta*w2;
            flag_tmp(j) = 1;
        else
            c_tmp(j)    = -1.0e+25;
            flag_tmp(j) = 0;
        end
        j = j+1;
    end
end
    
[c pos] = max(c_tmp);

w(1) = w_tmp(pos,1);
w(2) = w_tmp(pos,2);

flag = flag_tmp(pos);
=======
function [w,c,flag]=inner_fun(h,nv,v,pi,w_min,pi_star,delta)

% Choosing future values for each action
% v are the vertices (reduced before)
% nv number of vertices
% h search subgradient
% current payoff

alpha = 0:.01:1;
na = length(alpha);
T = (nv-1)*na;

c_tmp = zeros(T,1);
w_tmp = zeros(T,2);
flag_tmp = zeros(T,1);

j=1;
for aa=1:na
    for ih=1:nv-1
        w1 = (1-alpha(aa))*v(ih,1)+alpha(aa)*v(ih+1,1);
        w2 = (1-alpha(aa))*v(ih,2)+alpha(aa)*v(ih+1,2);
        
        ic_1 = (1-delta)*pi(1)+delta*w1-(1-delta)*pi_star(1)-delta*w_min(1);
        ic_2 = (1-delta)*pi(2)+delta*w2-(1-delta)*pi_star(2)-delta*w_min(2);
        if ic_1>=0 & ic_2>=0
            c_tmp(j)    = h(1) * ((1-delta)*pi(1)+delta*w1) + h(2) * ((1-delta)*pi(2)+delta*w2);
            w_tmp(j,1)  = (1-delta)*pi(1)+delta*w1;
            w_tmp(j,2)  = (1-delta)*pi(2)+delta*w2;
            flag_tmp(j) = 1;
        else
            c_tmp(j)    = -1.0e+25;
            flag_tmp(j) = 0;
        end
        j = j+1;
    end
end
    
[c pos] = max(c_tmp);

w(1) = w_tmp(pos,1);
w(2) = w_tmp(pos,2);

flag = flag_tmp(pos);
>>>>>>> 0b5165c8e0cbd0e4f1dde7640d267055799a76ba
