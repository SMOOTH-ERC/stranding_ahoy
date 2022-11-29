
%%   Transition expectations model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the model file: 
% - Contains only Probit model features in s and all other variables are in t
% - It is defined as a function, with key parameters being specified in the run file
% - The function returns an 'output' containing all the simulation results
% - The output is then loaded back into the run file, where results are analysed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[output] = Analytical_function_EC_15Jan21(T,S, r_log_sigma, SigmaMax, SigmaMin, networkexternality, rho,  theta, speed, initial_L, initial_H, g_xi, theta_xi, g_FF, alpha, betad) %%betad relevant for Stable Distributions                                  

%% Preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latest version for the Campiglio et al. (2022) Probit model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was created by LD on 11/21/2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preallocation of variables

% Electricity
e_d             = NaN(T+1+1,1);       % Total electricity demand 
e_ed            = NaN(T+1+1,T+1+S,1); % Total electricity demand 
e_L             = NaN(T+1,1);         % Electricity output by low-carbon technology 
e_s_L           = NaN(T+1,1);         % Electricity supply capacity by low-carbon technology
e_L_e           = NaN(T+1,T+1+S,1);   % Expected electricity output by low-carbon technology 
e_H             = NaN(T+1,1);         % Electricity output by high-carbon technology 
e_s_H           = NaN(T+1,1);         % Electricity supply capacity by high-carbon technology 
e_H_e           = NaN(T+1,T+1+S,1);   % Expected electricity output by high-carbon technology 
p_E             = NaN(T+1,1);         % Electricity market price


% Rates of low-carbon electricity related variables in totals
ell_E           = NaN(T+1,1);         % Share of low carbon electricity production in total prod.
ell_E_e         = NaN(T+1,T+1+S,1);   % Expected share of low carbon electricity production in total production
ell_K           = NaN(T+1,1);         % Share of low carbon electricity capital in total cap.
ell_I           = NaN(T+1,1);         % Share of low carbon electricity investment in total inv

% Capital stock and investment
K_L             = NaN(T+1,1);         % Low-carbon capital stock 
K_H             = NaN(T+1,1);         % High-carbon capital stock
K_H_d           = NaN(T+1,T+1+S,1);   % Desired high-carbon capital stock in the future
K_d             = NaN(T+1,1);         % Desired total capital stock
i_d             = NaN(T+1,1);         % Total Investment demand
i_L             = NaN(T+1,1);         % Investment in low-carbon capital
i_H             = NaN(T+1,1);         % Investment in high-carbon capital

dep_L           = NaN(T+1,1);         % Depreciation of Low-carbon capital vintaged stock 
dep_H           = NaN(T+1,1);         % Depreciation of high-carbon capital vintaged stock 
xi_L            = NaN(T+1,1);         % Productivity low-carbon capital
xi_H            = NaN(T+1,1);         % Productivity high-carbon capital
u_L             = NaN(T+1,1);         % Low-carbon capacity utilization rate 
u_H             = NaN(T+1,1);         % High-carbon capacity utilization rate 
ck_L            = NaN(T+1,1);         % Low-carbon unit capital costs 
ck_H            = NaN(T+1,1);         % High-carbon unit capital costs
ck_L_e          = NaN(T+1, T+1+S, 1); % Expected Low-carbon unit capital costs 
ck_H_e          = NaN(T+1, T+1+S, 1); % Expected High-carbon unit capital costs

% Fossil fuels
ff              = NaN(T+1,1);           % Fossil fuel input to high-carbon electricity sector (real)
p_FF            = NaN(T+1,1);           % Price of fossil fuels
p_FF_e          = NaN(S+1,1);           % Expected price of fossil fuels

% Loans 
L_L             = NaN(T+1,1);           % Debt high-carbon electricity sector
L_H             = NaN(T+1,1);           % Debt low-carbon electricity sector
rf_int          = NaN(T+1,1);           % Risk-free interest rate (target interest rate set by the central bank)
L_markup_L      = NaN(T+1,1);           % Mark-up rate on loans for low carbon electricity 
L_markup_H      = NaN(T+1,1);           % Mark-up rate on loans for high carbon electricity
int_L           = NaN(T+1,1);           % Interest rate on loans to low-carbon electricity sector
int_H           = NaN(T+1,1);           % Interest rate on loans to high-carbon electricity sector
alpha_L         = NaN(T+1,1);           % Cost recovery factor low-carbon capital
alpha_H         = NaN(T+1,1);           % Cost recovery factor high-carbon capital

% Expected utilisation
mu_u_H          = NaN(T+1,T+1+S,1);     % Central expectations of future high-carbon utilisation rate
sigma_u_H       = NaN(T+1,T+1+S,1);     % Uncertainty on utilisation rate expectations (standard deviation of error term)
sigma2_uH       = NaN(T+1,T+1+S,1);     % Variance of utilisation error term in time s

% Expected profit rates
mu_pi_H         = NaN(T+1,T+1+S,1);     % Central expectations of future high-carbon profit rate
mu_pi_L         = NaN(T+1,T+1+S,1);     % Central expectations of future low-carbon profit rate
mu_pi_H_disc    = NaN(T+1,T+1+S,1);     % Discounted expected high-carbon profit rates
mu_pi_L_disc    = NaN(T+1,T+1+S,1);     % Discounted expected low-carbon profit rates
sigma2_pi_H     = NaN(T+1,T+1+S,1);     % Variance of profit rate error term in time s
gamma_H         = NaN(T+1,1);           % Synthetic variable =xi_H(t)*(p_E(t) - p_FF(t)/xi_F)
gamma_H_e       = NaN(T+1,T+1+S,1);     % Expected synthetic variable
dr_H            = NaN(T+1,1);           % Debt repayment on a high-carbon capital unit
dr_L            = NaN(T+1,1);           % Debt repayment on a low-carbon capital unit
lower_bound_pi_h= NaN(T+1,T+1+S,1);     % Lower bound for high carbon profit rate 
upper_bound_pi_h= NaN(T+1,T+1+S,1);     % Upper Bound for low-carbon profit rate

%Actual profit Rates
Pi_L            = NaN(T+1,1);           % Actual low-carbon profit rate
Pi_H            = NaN(T+1,1);           % Actual high-carbon profit rate

% Expected returns
mu_r_L          = NaN(T+1,1);           % Central expectation of low-carbon return rate
mu_r_H          = NaN(T+1,1);           % Central expectation of high-carbon return rate
sigma2_r_H      = NaN(T+1,1);           % Variance of return rate error term
sigma_r_H       = NaN(T+1,1);           % Standard deviation of return rate error term
sigma_r_H_alt   = NaN(T+1,1); 
upper_bound_r   = NaN(T+1,1);           % Upper bound for high-carbon return rate
lower_bound_r   = NaN(T+1,1);           % Lower bound for high-carbon return rate
mu_phi          = NaN(T+1,1);           % Central expectations of return rate differential =mu_r_L(t)- mu_r_H(t)
phi             = NaN(T+1,1);           % mu_phi adjusted by sigma_r_H (argument of the CDF function)

% Annexes - Stable Distribution
beta_alt        = NaN(T+1,T+S+1,1);     % Skewness parameter of Stable Distribution 
beta_alt_dis    = NaN(T+1,T+S+1,1);     % Sum of skewness parameters with Stable Distribution
sigma2_pi_alt  = NaN(T+1,T+1+S,1);      % Variance - Profit rate - Stable ditribution
sigma_u_H_alt   =  NaN(T+1,1) ;         % Standard errors - Utilisation rate - Stable Distribution
ell_I_alt       = NaN(T+1,1);           % Share of low carbon electricity investment in total inv - Stable distribution

% Imagined Future Varables
K_H_e           = NaN(T+1, T+S+1);      % Expected High Carbon capital
xi_L_e          = NaN(S+1, 1);          % Expected renewable productivity

%% Parameters

% Demand growth
g_e             = 0.012;                 % Growth in electricity demand        

% Technical parameters       
xi_F            = 1/8.75;               % Productivity of fossil fuels
ck_L(:)         = 2.9;                  % Low-carbon capital unit cost
ck_H(:)         = 3.5;                  % High-carbon capital unit cost
u_H_tau         = 0.85;                 % Technical utilization rate

% Depreciation rates
delta_KL        = 0.04;                 % Depreciation rate low-carbon capital 
delta_KH        = 0.025;                % Depreciation rate high-carbon capital 

% Financial parameters: Interest rates and debt financing
psi_L           = 0.75;                 % Debt-to-investment ratio in low-carbon sector 
psi_H           = 0.7;                  % Debt-to-investment ratio in high-carbon sector 
LT              = 15;                   % Loan tenure rate

% Expectation parameters
beta            = 1/(1+rho);           % Discount factor
lb_uH           = 0;                   % Lower-bound for high-carbon capital utilization rate
hb_uH           = 1;                   % Higher-bound for high-carbon capital utilization rate

%Prices
p_FF(1)      = 0.0025;                 %Fossil fuel price billion $/ trillion Btu
p_E(1)       = 0.2;                    % billion$/TWh

%% Initial values 
% Utilisation 
u_L(1)          = 0.85;                % Initial utilisation rate of low-carbon capital (=1 due to priority feed-in of renewables)
u_H(1)          = 0.85;                % Initial utilisation rate of high-carbon capital (=0.85 to include spare capacity to keep the electricity system stable) 

% Capital stocks
K_L(1)        = initial_L;             % Initial low-carbon capital stock
K_H(1)        = initial_H;             % Initial high-carbon capital stock depends on the residual of total energy production (renewable priority feed-in) and the initial productivity and capacity utilization
K_d(1)        = K_L(1)+K_H(1);         % For simplicity, assume for first modelling period that desired total capital stock is equal to sum of actual capital stocks 

% Electricity Demand
e_d(1)        = 3243;                  % Initial energy demand (=3000 TWh) 
ell_E(1)      = 0.22;
xi_H(:)       = ((1-ell_E(1))*e_d(1) )/(u_H(1)*K_H(1));%6.132;         % Productivity of high-carbon capital
xi_L(:)       = (ell_E(1)*e_d(1) )/(u_L(1)*K_L(1));%2.628;             % Productivity of low-carbon capital


% Dynamic expected productivity - IRRELEVANT HERE
xi_L_e(1)       = xi_L(1);             % Initital Expected Productivty
xi_L_o          = xi_L_e(1)-0.01;           % Setting lower asymptote, here equal to productivity
for i=2:S
     xi_L_e(i) = xi_L_o+(xi_L_e(i-1)-xi_L_o)*(1+g_xi*(1-xi_L_e(i-1)/theta_xi)); %Increases logistically to theta_xi. g_xi is 0 in the current version, such that the producitivity is expected to be constant.
end

%Dynamic expected fossilfuel price - g_FF = 0 
p_FF_e(1)       = p_FF(1);             %
for i=2:S
     p_FF_e(i) = p_FF_e(i-1)*(1+g_FF); %
end

% Simplifying assumption: loans are initialized as debt-financed rate of capital stock
L_L(1)          = psi_L*K_L(1);         % Initial levels of low-carbon loans 
L_H(1)          = psi_H*K_H(1);         % Initial levels of high-carbon loans 
rf_int(:)       = 0.01;                 % risk-free interest rate (target interest rate set by the central bank )
L_markup_H(:)   = 0.035;                % Mark-up rate on loans for high-carbon electricity
L_markup_L(:)   = L_markup_H(:)*0.84;   % Mark-up rate on loans for low-carbon electricity 


%%  Equations

% BEGINNING TOTAL MODELLING LOOP
for t = 1:T  %Here T=1 for static version of the model
p_FF(t)      = 0.0025;                 %Fossil fuel price billion $/ trillion Btu
p_E(t)       = 0.2;                    % billion$/TWh


%% Electricity system 
  e_d(t+1)      =   e_d(t)*(1+g_e);  
  % Energy demand expectations are always correct
  for s=1:T+S-1
    e_ed(t,t)       =  e_d(t);  
    e_ed(t,s+1)     =  e_ed(t,s)*(1+g_e);  
  end

  % Electricity capacity evolution
  e_s_L(t)     = K_L(t)*xi_L(t);
  e_s_H(t)     = K_H(t)*xi_H(t);

  % Low-carbon electricity production
  if        K_L(t)== 0          
                e_L(t) = 0; 
  elseif    e_d(t)<= e_s_L(t)   
                e_L(t) = e_d(t) ;
  elseif    e_d(t)>  e_s_L(t)    
                e_L(t) = u_L(1)*e_s_L(t);
  end

  % High-carbon electricity production
  if        K_H(t)==0                          
                e_H(t) = 0;
  elseif    e_d(t)<=e_s_L(t)                    
                e_H(t) = 0;
  elseif    e_s_L(t)<e_d(t)<(e_s_L(t)+e_s_H(t))   
                e_H(t) = e_d(t)  - u_L(1)*e_s_L(t);
  elseif    e_d(t)>=e_s_L(t)+e_s_H(t)            
                e_H(t) = e_s_H(t);
  end
  % Share of low-carbon energy
  ell_E(t) = e_L(t)/(e_L(t) + e_H(t));

  % Capacity utilization rate
  u_L(t) = e_L(t)/e_s_L(t);
  u_H(t) = e_H(t)/e_s_H(t);
  % Fossil fuel consumption
  ff(t)  =  e_H(t)/xi_F;
  
%% Prices
  gamma_H(t)   = (p_E(t)-p_FF(t)/xi_F)*xi_H(t);                          % Synthetic variable to improve readability
  % Expected costs associated with change in fossil fuel price
  for s=1:S
  gamma_H_e(t,s) = (p_E(t)-p_FF_e(s)/xi_F)*xi_H(t);                      %Fossil fuel price constant in expectations in current version
  end
  %Current share of low-carbon capacity
  ell_K(t)      = K_L(t)/(K_L(t) + K_H(t));

  %Updating productivity with network externalities (default set to 0) - %IRRELEVANT HERE
  if networkexternality == 0
  ck_L(t+1)         = ck_L(1); 
  elseif networkexternality ==1
  ck_L(t+1)         = ck_L(1)*((1 + 0.5*ell_K(1)^2))/((1 + 0.5*ell_K(t)^2));
  end
  
%% Loans
  int_L(t)      = rf_int(t) + L_markup_L(t);                            % Interest rate on low-carbon loans is a mark-up on rf_int  
  int_H(t)      = rf_int(t) + L_markup_H(t);                            % Interest rate on high-carbon loans is a mark-up on rf_int
  alpha_L(t)    = (int_L(t)*(1+int_L(t))^LT )/((1+int_L(t))^LT - 1);    % Capital recovery factor for low-carbon loans
  alpha_H(t)    = (int_H(t)*(1+int_H(t))^LT )/((1+int_H(t))^LT - 1);    % Capital recovery factor for high-carbon loans
  % Debt costs
  dr_H(t)       = alpha_H(t)*psi_H*ck_H(t);                             % Expected high-carbon debt repayment in profitability assessment
  dr_L(t)       = alpha_L(t)*psi_L*ck_L(t);                             % Expected low-carbon debt repayment in profitability assessment
  % Profits 
  Pi_L(t)       = p_E(t)*xi_L(t)*u_L(t)-dr_L(t);                        %Expected total profit - Low-carbon
  Pi_H(t)       = (p_E(t)-(p_FF(t)/xi_F))*xi_H(t)*u_H(t)-dr_H(t);       %Expected total profit - High-carbon

%% Depreciation
dep_L(t) = delta_KL*K_L(t);
dep_H(t) = delta_KH*K_H(t);

%% Central stranding
ell_E_e(t,1:t)=ell_E(1:t);
e_L_e(t,1:t)=e_L(1:t);
e_H_e(t,1:t)=e_H(1:t);
K_H_e(t,1:t)=K_H(1:t);
for s=1:S
    ell_E_e(t,t+s)  =   ell_E_e(t,t+s-1)*(1+speed*(1-ell_E_e(t,t+s-1)/theta));                     % Expected low-carbon energy share
    e_L_e(t,t+s)    =   e_L_e(t,t+s-1)*(1+g_e)*(1+speed*(1-e_L_e(t,t+s-1)/(theta*e_ed(t,t+s-1)))); %Expected low-carbon energy demand
    e_H_e(t,t+s)    =   e_ed(t,t+s)-e_L_e(t,t+s);                                                  %Expected high-carbon energy demand
    K_H_d(t,t+s)    =   e_H_e(t,t+s)/(u_H(1)*xi_H(t));                                             %Expected demand for high-carbon capital in case of no stranding
    K_H_e(t,t+s)    =   max(K_H_d(t,t+s),(1-delta_KH)*K_H_e(t,t+s-1));                             %Actual expected demand for high-carbon capital
    mu_u_H(t,t+s)   =   e_H_e(t,t+s)/(xi_H(t)*K_H_e(t,t+s));                                       %Expected utilisation rate - Central Expectation
end



%% Uncertainty of Expectations
  for s = 1:S
    % Uncertainty around mu_u_H
    sigma_u_H(t,t)          = SigmaMin;
    if (SigmaMax>=SigmaMin)
    sigma_u_H(t,t+s)        = (sigma_u_H(t,t+s-1) +  r_log_sigma*sigma_u_H(t,t+s-1)*(1-sigma_u_H(t,t+s-1)/(SigmaMax)));
    else
    sigma_u_H(t,t+s)        = 0;    
    end
    sigma2_uH(t,t+s-1)      = sigma_u_H(t,t+s-1)^2;
    % Expected unit profits mu_pi
    mu_pi_H(t,t+s)          =  gamma_H_e(t,s)*mu_u_H(t,t+s)- dr_H(t);
    mu_pi_H_disc(t,t+s)     =  beta^s*mu_pi_H(t,t+s);
    mu_pi_L(t,t+s)          =  p_E(t)*xi_L_e(s)*u_L(t) - dr_L(t);
    mu_pi_L_disc(t,t+s)     =  beta^s*mu_pi_L(t,t+s);
    % Uncertainty around expected mu_pi_H
    sigma2_pi_H(t,t)    = sigma_u_H(t,t)*(gamma_H(t))^2;
    sigma2_pi_H(t,t+s)    = (beta^(s)*sigma_u_H(t,t+s)*(gamma_H(t)))^2;
    sigma2_pi_alt(t,t)    = (sigma_u_H(t,t)*(gamma_H(t)))^alpha;
    sigma2_pi_alt(t,t+s)    = (beta^(s)*sigma_u_H(t,t+s)*(gamma_H(t)))^alpha;
    % Lower and upper bounds for expected unit profits
    lower_bound_pi_h(t,t+s) =  beta^s*(gamma_H(t)*lb_uH - dr_H(t));
    upper_bound_pi_h(t,t+s) =  beta^s*(gamma_H(t)*hb_uH - dr_H(t));
  end
   
  
%% Expected return rates
  % Expected return rates
  mu_r_L(t)       = sum(mu_pi_L_disc(t,t+1:t+S));
  mu_r_H(t)       = sum(mu_pi_H_disc(t,t+1:t+S));
  % mu_phi is the difference between the two expected return rates
  mu_phi(t)       = mu_r_L(t)- mu_r_H(t);
  % Uncertainty around mu_r_H
  sigma2_r_H(t)   = sum(sigma2_pi_H(t,t+1:t+S));
  sigma_u_H_alt(t) = sum(sigma2_pi_alt(t,t+1:t+S));
  % Defining the shifting/rescaling of the underlying normal distribution for the Probit 
  sigma_r_H(t)    = sqrt(sigma2_r_H(t));
  sigma_r_H_alt(t)  = sigma_u_H_alt(t)^(1/alpha);
  % Defining bounds for profits (Note that the names are reversed compared to those of u because they are subtracted from R_L
  upper_bound_r(t)= (mu_r_L(t) - sum(lower_bound_pi_h(t,t+1:t+S)))/sigma_r_H(t); 
  lower_bound_r(t)= (mu_r_L(t) - sum(upper_bound_pi_h(t,t+1:t+S)))/sigma_r_H(t);
  % Synthetic Variable: argument of the CDF
  phi(t)             = mu_phi(t)/(sigma_r_H(t));


  %%%Adding Investment rule for Stable Distribution - Neutralised in vanilla version
  %Defining dynamic beta (asymptote) for Stable Distribution
 %   beta_alt(t,t) = -1;
 % for s= 1:S
 %     beta_alt(t,t+s) = 1-mu_u_H(t,t+s)*2/u_H_tau;
 %  end
%%%%%Defining skewness parameter  
%   beta_alt_dis(t,t) = betad;
%  beta_alt_dis(t,t+1) = beta_alt(t,t+1);
% for s=2:S
%    beta_alt_dis(t,t+s) = (beta_alt(t,t+s-2)*sigma2_pi_alt(t,t+s-2) + beta_alt_dis(t,t+s-1)*sigma2_pi_alt(t,t+s-1))/(sum(sigma2_pi_alt(t,t:t+s)));
% end
%pd = makedist('Stable','alpha',alpha,'beta',beta_alt_dis(t,t+S),'gam',sigma_r_H_alt(t),'delta',0)  ;%Definining stable distribution density based on parameters
%ell_I_alt(t) = ((1-cdf(pd,-mu_phi(t)/gamma_H(t))) - (phi(t) > lower_bound_r(t))*(phi(t) < upper_bound_r(t))*(1-cdf(pd2,-lower_bound_r(t)/gamma(t))) - (phi(t) > upper_bound_r(t))*(1-cdf(pd,-upper_bound_r(t)/gamma(t))))/((1-cdf(pd,-upper_bound_r(t)/gamma(t))) - (1-cdf(pd,-lower_bound_r(t)/gamma(t))));


%% INVESTMENT DECISION     
  %Probability to shift from low- to high-carbon investment
  if phi(t) == -Inf %Indicating asymptotic behaviour
  ell_I(t) = 0;
  elseif phi(t) < upper_bound_r(t)
  ell_I(t) = (cdf('Normal',phi(t),0,1) - (phi(t) > lower_bound_r(t))*(phi(t) < upper_bound_r(t))*cdf('Normal',lower_bound_r(t),0,1)  -  (phi(t) > upper_bound_r(t))*cdf('Normal',upper_bound_r(t),0,1))/(cdf('Normal',upper_bound_r(t),0,1) - cdf('Normal',lower_bound_r(t),0,1)) ;
  elseif phi(t) > upper_bound_r(t) || phi(t) == Inf %Indicating asymptotic  and out-of-range behaviour
  ell_I(t)=1;
  end
 
%% TOTAL DESIRED INVESTMENT - Neutralised since we do not use chronological time
%   xi_e(t)     = u_L(t)*xi_L(t)*ell_K(t)+ mu_u_H(t,t+1)*xi_H(t)*(1-ell_K(t)) ; 
%  %u_bar(t)   = u_L(t)*ell_K(t) +u_H(t)*(1-ell_K(t)) ;
%   Gamma1      = (1-ell_I(t))*K_L(t) + ell_I(t)*(dep_L(t) + dep_H(t) - K_H(t)) - dep_L(t);
%   Gamma2      = ell_I(t)*K_H(t) + (1-ell_I(t))*(dep_L(t) + dep_H(t) - K_L(t)) - dep_H(t);
%   % Desired  capital stock
%   K_d(t+1)    = (e_d(t+1) - xi_L(t)*Gamma1 - u_H_tau*xi_H(t)*Gamma2)/(u_H_tau*xi_H(t)*(1-ell_I(t)) + xi_L(t)*ell_I(t));
%  
% %   % Desired  capital stock
% %     K_d(t+1ell_I)    = e_d(t+1)/(xi_L(t)*(1-h_e(t,t+1))+xi_H(t)*u_H_tau*h_e(t,t+1));
%   
%   % Desired investments
%   i_d(t)      = max(K_d(t+1)-(K_L(t)+K_H(t)),0) + dep_L(t) + dep_H(t);
%   
%   % Actual investment in technologies - where it all comes together  : PROBIT model
%   i_L(t)    =  i_d(t).*ell_I(t); 
%   i_H(t)    =  i_d(t) - i_L(t);

  %% ACCOUNTING
  % Capital stock laws of motion
  K_L(t+1)      =   i_L(t)+ K_L(t) - dep_L(t);
  K_H(t+1)      =   i_H(t)+ K_H(t) - dep_H(t);
  % Loan evolution
  L_L(t+1)      = L_L(t) - alpha_L(t)*L_L(t) + psi_L * i_L(t);
  L_H(t+1)      = L_H(t) - alpha_H(t)*L_H(t) + psi_L * i_H(t);
  
     
% END TOTAL MODELLING LOOP  

end


%% Results
res_t=table([1:t+1]',ell_I,ell_K,ell_E, K_d, u_H, phi, mu_phi, sigma_r_H, mu_r_L, mu_r_H, dr_H,dr_L);
res_t.Properties.VariableNames(1)={'time_t'};

t_s=1;
res_s=table([t_s:t_s+s]',mu_u_H(t_s,t_s:t_s+S)',sigma_u_H(t_s,t_s:t_s+S)',mu_pi_H(t_s,t_s:t_s+S)',mu_pi_L(t_s,t_s:t_s+S)');
res_s.Properties.VariableNames = {'time_s','mu_u_H' 'sigma_u_H', 'mu_pi_H', 'mu_pi_L'};

save('Results')
output='Results';
end

