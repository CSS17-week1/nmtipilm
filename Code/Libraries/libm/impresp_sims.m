% Computes impulse responses for solution of linear RE models,
% assumes i.i.d. shock process (z)
% Inputs:
%   G1 and Impact: Output from Sims' package
%   IMP_SELECT:    indx of variables for which to draw impulse responses
%   HORIZON:       number of periods for which to compute IR
%   PERIOD:        number of model periods per year
%   varnames:      names of variables
%   shocknames:    names of shocks
% Ouput:
%   Cell array of impulse responses, one element for each shock
% Plots impulse responses for the variables IMP_SELECT, each shock in one figure
% Adapted from impulse-response routine in Uhlig's toolkit
% Michael Reiter, October 2005 
% (current version very messy)
function Resp_mat = impresp_sims(G1,Impact,IMP_SELECT,HORIZON,PERIOD,varnames,shocknames,scale_shock,scale_y);
  global file_suffix;
  global irsims_doprint;

  if(~exist('irsims_doprint'))
    irsims_doprint = 1;
  end;



  nz = size(Impact,2);
  if(nargin<9 | isempty(scale_y))
    scale_y = ones(size(G1,1),1);
  end;
  if(nargin<8 | isempty(scale_shock))
    scale_shock = ones(nz,1);
  end;
  NN = zeros(nz,nz);

  [m_states,k_exog] = size(Impact);
  n_endog = 0;
  Time_axis = (0 : (HORIZON-1))/PERIOD;

  %  Resp_mat = zeros(m_states+n_endog+k_exog,HORIZON,k_exog);
  t_reverse = 0;
  Resp_mat = cell(k_exog,1);
  for shock_counter = 1 : k_exog;
    Response = zeros(m_states+n_endog+k_exog,HORIZON);
    iVar = m_states+n_endog+shock_counter;
    Response(iVar,1) = 1;
    II_lag = [ G1, zeros(m_states,n_endog),zeros(m_states,k_exog)
%              RR, zeros(n_endog, n_endog),zeros(n_endog, k_exog)
              zeros(k_exog,(m_states+n_endog)), NN                ];
    aux1 = [ zeros(m_states,(m_states+n_endog)), Impact];
%    aux2 = [ zeros(n_endog, (m_states+n_endog)), SS];
    aux3 = [ zeros(k_exog,  (m_states+n_endog)), zeros(k_exog,k_exog) ];
%    II_contemp = eye(m_states+n_endog+k_exog) + [aux1;aux2;aux3];
    II_contemp = eye(m_states+n_endog+k_exog) + [aux1;aux3];
    % describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';
    Response(:,1) = II_contemp*Response(:,1);
    for time_counter = 2 : HORIZON;
      Response(:,time_counter) = II_contemp*II_lag*Response(:,time_counter-1);
      if(time_counter==t_reverse)
	Response(iVar,t_reverse) = Response(iVar,t_reverse) - 1;
      end;
    end;
    Response = Response(1:m_states,:);
    Response = Response*scale_shock(shock_counter);
    %assert(length(scale_y)==size(Response,1));
    indx_scale = scale_y ~= 0;
    Response = Response .* repmat(scale_y,1,size(Response,2));
    Resp_mat{shock_counter} =  Response;
    fname = sprintf('IRshock%d%s.txt',shock_counter,file_suffix);
%    disp(sprintf('Impulse response to %s   in file %s',varnames(iVar,:),fname));
%    save_ascii(fname,varnames,Resp_mat');
%    save_ascii(fname,Resp_mat');

  if(irsims_doprint)
    figure(shock_counter);
    tt = (0:HORIZON-1) / PERIOD;
    plot(tt,Response(IMP_SELECT,:)');
    legend(varnames(IMP_SELECT,:));
    title(['impulse response to ' shocknames(shock_counter,:)]);
    xlabel('years after impulse');
    ylabel('deviations from steady state');
  end
  end;

  return;
  for state_counter = 1 : m_states
    Response = zeros(m_states+n_endog+k_exog,HORIZON);
    Response(state_counter,1) = 1;
    II_lag = [ G1, zeros(m_states,n_endog),zeros(m_states,k_exog)
%              RR, zeros(n_endog, n_endog),zeros(n_endog, k_exog)
              zeros(k_exog,(m_states+n_endog)), NN                ];
    II_contemp = eye(m_states+n_endog+k_exog) + ...
        [ zeros(m_states,(m_states+n_endog)), Impact
%         zeros(n_endog, (m_states+n_endog)), SS
         zeros(k_exog,  (m_states+n_endog)), zeros(k_exog,k_exog) ];
    % describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';
    %    Response(:,1) = II_contemp*II_lag*Response(:,1); %MR: don't know why he is doing that!!
    for time_counter = 2 : HORIZON;
      Response(:,time_counter) = II_contemp*II_lag*Response(:,time_counter-1);
    end;
    fname = sprintf('IRstate%d%s.txt',state_counter,file_suffix);
%    disp(sprintf('Response to change in %s in file %s',varnames(state_counter,:),fname));
%    save_ascii(fname,varnames,Response');
  end;
