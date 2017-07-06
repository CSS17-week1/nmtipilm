
%==========================================================================     
% Computing Set of SubGame Perfect Equilibria                             
% Infinitely Repeated Prisoner's Dilemma                                  
% Based on Computing Supergame Equilibria Conklin, Judd and Yeltekin (2003)
%
% Author: Pablo D'Erasmo
% Date: May 2007
%==========================================================================

clear all
tic

% Number of subgradients
D = input('Number of search subgradients (> 3 ) = ');            

% Discunt Factor
delta = 0.8;

%--------------------------------------------------------
% Matrix of stage game payoffs (columns are players)
% moves clockwise
%--------------------------------------------------------

sgp = zeros(4,2);
sgp(1,1) = 4; sgp(1,2) = sgp(1,1);   
sgp(3,1) = 2; sgp(3,2) = sgp(3,1);
sgp(2,1) = 0; sgp(4,2) = sgp(2,1);
sgp(4,1) = 6; sgp(2,2) = sgp(4,1);

fprintf('\nStage Game Payoffs:\n\n\t\t\tCOOP.\tDEFECT\n')
fprintf('\tCOOP.  (%g,%g)\t (%g,%g)\n',sgp(1,1),sgp(1,2),sgp(2,1),sgp(2,2))
fprintf('\tDEFECT (%g,%g)\t (%g,%g)\n\n',sgp(4,1),sgp(4,2),sgp(3,1),sgp(3,2))

%----------------------------------------------------------------------
% Stage game deviations (columns are players)
% rows are initial position in the payoff matrix
% Some of them are trivial but just for completion
%----------------------------------------------------------------------

pi_star(1,1) = sgp(4,1);
pi_star(2,1) = sgp(3,1);
pi_star(3,1) = sgp(3,1);
pi_star(4,1) = sgp(4,1);

pi_star(1,2) = sgp(2,2);
pi_star(2,2) = sgp(2,2);
pi_star(3,2) = sgp(3,2);
pi_star(4,2) = sgp(3,2);

% -------------------------------------------------------
% Generating vector of directions 
% -------------------------------------------------------

jump = 360/D;
theta = 0:jump:360;
% Vector of subgradients is H 
H = [(cosd(theta))' (sind(theta))'];

% -------------------------------------------------------
% Initial bounday points
% -------------------------------------------------------
% Creating middle points and expanding
r = 8;
z = r.*(H(1:D,:)+H(2:D+1,:))/2;

% Last vector of H is redundant
H = H(1:D,:);

% -------------------------------------------------------
% Initial Vector of levels
% -------------------------------------------------------
c = sum(H.*z,2);

%---------------------------------------------------------
% Approximating the set of equilibrium Payoffs 
%---------------------------------------------------------

% Initial Points
wo=z;
options = optimset('fmincon');
options = optimset(options,'LargeScale','on','Algorithm','sqp','Display','off','MaxIter',10000);
tol_w = 1.0e-06;

%---------------------------------------------------------
% Outer approximation
%---------------------------------------------------------
disp('=====================================================================')

distw=1;
HH = [H;-delta 0;0 -delta];  % Augmenting the set of constraints to include IC 
while distw>tol_w
    % Computing the min payments from initial set
    w_min(1) = min(wo(:,1));
    w_min(2) = min(wo(:,2));
    
    for ih=1:D % for all search directions
        for ai=1:4 % for all possible strategy profiles
            cc = [c;-((1-delta)*pi_star(ai,1)+delta*w_min(1)-(1-delta)*sgp(ai,1)); -((1-delta)*pi_star(ai,2)+delta*w_min(2)-(1-delta)*sgp(ai,2))]; % Augmenting the set of constraints to include IC 
            x0=[2;2];
            [w_tmp(ai,:),c_tmp(ai),flag]   = fmincon('outer_fun',x0,HH,cc,[],[],[],[],[],options,H(ih,:),delta,sgp(ai,:));
        
            if flag~=1      % No w satisfying constraints
                c_tmp(ai)=-1e+35;
            else
                c_tmp(ai)=-c_tmp(ai);
            end
        end
    
        [c_new(ih), pos] = max(c_tmp); % get the largest c and the index for its position
        wo_new(ih,:) = (1-delta)*sgp(pos,:)+delta*w_tmp(pos,:); % calculate the associated payoff
    end
    
    wo = wo_new;
    distw = max(abs(c_new'-c));
    disp(['Outer Approximation Distance = ', num2str(distw)])
    c = c_new';    
end
            
disp('=====================================================================')


%---------------------------------------------------------
% Inner approximation
%---------------------------------------------------------

% Shrinking the outer approximation set by 3 %
for ih=1:D % for all search directions
    if wo(ih,1)>=min(wo(:,1))-1.e-08 & wo(ih,1)<=min(wo(:,1))+1.e-08 
        wi(ih,1) = wo(ih,1);
    else
        wi(ih,1) = wo(ih,1)-.02*abs(wo(ih,1));
    end
    if wo(ih,2)>=min(wo(:,2))-1.e-08 & wo(ih,2)<=min(wo(:,2))+1.e-08
        wi(ih,2) = wo(ih,2);
    else
        wi(ih,2) = wo(ih,2)-.02*abs(wo(ih,2));
    end
end



distw=1;
iteri = 1;
while distw>tol_w
    % Computing the min payments from initial set
    w_min(1) = min(wi(:,1));
    w_min(2) = min(wi(:,2));
    
    k = convhull(wi(:,1),wi(:,2));
    v = wi(k,:);
    nv = length(v);
    
    for ih=1:D % for each search direction
        for ai=1:4 % for each strategy profile
            [w_tmp(ai,:),c_tmp(ai),flag]   = inner_fun(H(ih,:),nv,v,sgp(ai,:),w_min,pi_star(ai,:),delta);
        
            if flag~=1      % No w satisfying constraints
                c_tmp(ai)=-1e+35;
            end
        end
    
        [c_new(ih), pos] = max(c_tmp);
        wi_new(ih,:)     = w_tmp(pos,:);
    end
    distw = max(max((wi-wi_new).^2));
    disp(['Inner Approximation Distance = ', num2str(distw)])
    wi = wi_new;
    iteri = iteri + 1;
end

disp('=====================================================================')

toc

% making the plots

na = 6;
for ih=1:D
    if H(ih,2)==0
        aa(ih,:) = ones(1,na)*wo(ih,1);
        bb(ih,:) = 1:na;
    else
        aa(ih,:) = 1:na;     
        bb(ih,:) = c(ih)/H(ih,2)-H(ih,1)/H(ih,2).*aa(ih,:);
    end
end

% outer approximation
figure;
for ih=1:D
    plot(wo(ih,1),wo(ih,2),'o',aa(ih,:),bb(ih,:),'-')
    title(['Outer approximation L = ', num2str(D)])
    axis([0.5 6.5 0.5 6.5])
    %pause
    hold on
end


% Outer Approximation vertices
vert=[];
for i=1:D-1
    aux=inv([H(i,:);H(i+1,:)])*[c(i);c(i+1)];
    vert=[vert;aux'];    
end
aux=inv([H(D,:);H(1,:)])*[c(D);c(1)];
vert=[vert;aux'];        

figure;
scatter(vert(:,1),vert(:,2));hold    
vert=[vert;vert(1,:)];     
plot(vert(:,1),vert(:,2));
title(['Outer approximation vertices L = ', num2str(D)])
axis([0.5 6.5 0.5 6.5]);
print -depsc2 -r300 f_0124;

% inner approximation
figure;
plot(v(:,1),v(:,2),'o-')
title(['Inner approximation L = ', num2str(D)])
axis([0.5 6.5 0.5 6.5])
print -depsc2 -r300 f_0224;

% inner and outer approximation
figure;
plot(v(:,1),v(:,2),'b.-',vert(:,1),vert(:,2),'ro-');
title(['Inner and Outer approximation L = ', num2str(D)])
legend('Inner', 'Outer')
axis([0.5 6.5 0.5 6.5])
print -depsc2 -r300 f_0324;
     



