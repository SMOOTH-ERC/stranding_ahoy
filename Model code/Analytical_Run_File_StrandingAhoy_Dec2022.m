
%%  Transition expectations model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the run file (Version LD-11/21/2022):
% - It defines the values of the most relevant scenario parameters
% - Parameter ranges for sensitivity analysis can be defined
% - It calls a function with all the model equations
% - It stores all results into a structure called sens_res
% - It calls on various scripts to plot figures (comment/uncomment as desired)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
 
%% Setting sensitivity analysis
% Specify up to two parameters for sensitivity analysis  
% Insert the parameter names as strings 
% (NB: the chosen parameters need to be redefined in the 'Setting parameters' section)
sens_par_name_1 = "Long-run belief heterogeneity ($\bar{\sigma}_{\pi_{H}}$)";
sens_par_name_2 = "\begin{tabular}{c} Short-run belief \\ heterogeneity ($\sigma_{u_{H},0}$)\end{tabular}" ;
sens_var_name   = "\begin{tabular}{c} Share of \\ low-carbon \\ investment \\ $(\ell_{I})$\end{tabular}";
% To copy when needed "Length of the planning horizon $S$", "Intrinsic growth rate$b_{\ell}$", "$Short-run belief heterogeneity \sigma_{0}$", "$Long-rin belief heterogeneity \bar{\sigma}$", "$Corporate discount rate \rho$", "$Carrying capacity \bar{\ell}$"
% Define the range of the sensitivity analysis 
sens_par_range_1 = 0.5:0.05:2;
sens_par_range_2 = 0:0.05:0.5;

%% Start of the loops  
% Initialize the second scenario count parameter
sens_run_2=  0;
% Start of the second sensitivity analysis loop
for sens_par_2 = sens_par_range_2    
    sens_run_2 = sens_run_2+1;
    % Initialize the first scenario count parameter
    sens_run_1 = 0;
    % Start of the first sensitivity analysis loop
    for sens_par_1= sens_par_range_1
        sens_run_1 = sens_run_1+1;  
         
%% Parameter setting
        %Ambitious/Unambitious switch for central expectation
        ambitious = 0;
        extended_uncertainty=1;
        % NB. If sensitivity analysis is run, the chosen parameters must be assigned the ranges specified above
        T                  = 1              ; % Total length of simulations 
        %% Sensitivity parameters
        S                  = 20                                   ; % Length of planning horizon        (default=20)
        speed              = ambitious*0.25 + (1-ambitious)*0.15   ; % Expected speed of logistic transition. To replace by sens_par_1 or sens_par_2 to reproduce Fig.5
        r_log_sigma        = 0.28                                  ; % Intrinsic growth rate of sigma_u_H - Fixed here for simplicity
        SigmaMin           = sens_par_2                           ; % Starting value for sigma_u_H      (default=0.01)
        SigmaMin_piL       = extended_uncertainty*0.01                                 ;  %Starting value for sigma_pi_L      (default=0 (No belief dispersion))
        SigmaMin_piH       = extended_uncertainty*0.01                                 ;  %Starting value for sigma_pi_H      (default=0 (No belief dispersion))
        SigmaMax           = sens_par_1                           ; % Carrying capacity of sigma_u_H    (default=1)
        SigmaMax_piL       = 0.5                                  ;  %Carrying capacity for sigma_pi_L      (default=1 )
        SigmaMax_piH       = 0.5                                 ;  %Carrying capacity for sigma_pi_H      (default=1 )
        rho                = 0.05                                  ; % Discount Rate (default = 0.05)
        theta              = 0.85                                   ; % Carrying Capacity of expected low carbon energy prod. share (default=0.85)
        
        %%%Stable distribution Annex
        alpha              = 2                                     ; %Kurtosis parameter for stable distribution (default = 2)
        betad              = 0                                     ; %Initial skewness parameter for stable distribution (default = 0)    

        %Further model features - IRRELEVANT IN CURRENT VERSION
        networkexternality = 0                                     ; % Switch for network externalities function (default = 0)
        theta_xi           = 3.8                                   ; % Expected Limit of Low-carbon Productivity  (default = 2.893?)
        g_xi               = 0                                     ; % Expected Speed of Productivity Growth (default = 0)
        g_FF               = 0                                     ; % Expected growth of fossil fuel price (Default = 0)

        %Initial values
        initial_L          = 288           ; % Initial Low carbon capital stock 
        initial_H          = 947 - initial_L           ; % Initial High carbon capital stock
        ell_E_start        = 0.22           ; % Initial low-carbon energy production
        

        %% Run model
        % Model equations are called as a function file
        output = Analytical_function_StrandingAhoy_Dec2022(T,S, r_log_sigma, SigmaMax, SigmaMin, networkexternality, rho, theta, speed, initial_L, initial_H, g_xi, theta_xi, g_FF, alpha, betad, SigmaMin_piL, SigmaMin_piH, SigmaMax_piL, SigmaMax_piH,extended_uncertainty, ell_E_start);
        load(output)
       
        %% Store results
        % Results are stored into a structure with specific fields for most relevant variables
        sens_res(sens_run_1,sens_run_2).name            =   "SensitivityRun_for_"+sens_par_name_1+"="+sens_par_1+"_and_"+sens_par_name_2+"="+sens_par_2+"";
        sens_res(sens_run_1,sens_run_2).ell_I           =   ell_I;
        sens_res(sens_run_1,sens_run_2).ell_I_alt       =   ell_I_alt;
        %sens_res(sens_run_1,sens_run_2).ell_I_alt2      =  ell_I_alt2;
        sens_res(sens_run_1,sens_run_2).ell_E           =   ell_E;
        sens_res(sens_run_1,sens_run_2).ell_E_e         =   ell_E_e;
        sens_res(sens_run_1,sens_run_2).ell_K           =   ell_K;
        sens_res(sens_run_1,sens_run_2).u_H             =   u_H;
        %sens_res(sens_run_1,sens_run_2).h_e             =   h_e;
        sens_res(sens_run_1,sens_run_2).mu_u_H          =   mu_u_H;
        sens_res(sens_run_1,sens_run_2).sigma_u_H       =   sigma_u_H;
        sens_res(sens_run_1,sens_run_2).mu_pi_H         =   mu_pi_H;
        sens_res(sens_run_1,sens_run_2).mu_r_L          =   mu_r_L;
        sens_res(sens_run_1,sens_run_2).mu_r_H          =   mu_r_H;
        sens_res(sens_run_1,sens_run_2).mu_phi          =   mu_phi;
        sens_res(sens_run_1,sens_run_2).phi             =   phi;
        %sens_res(sens_run_1,sens_run_2).change_E        =   change_E;
        %sens_res(sens_run_1,sens_run_2).change_K        =   change_K;
        sens_res(sens_run_1,sens_run_2).sigma_r_H       =   sigma_r_H;
% %         
%         plot([1:S+1], transpose(h_e))
%         hold on
% End of loops
    end
end


%%  Figures - Default is surface plot
var_sens_s = 'ell_I';
var_sens_s_latex = '$\ell_{I}$';
t_focus=1;
run('Figure_codes/Sensitivity_figures_in_time_s');

%%  Figures - Contour
var_sens_s = 'ell_I';
var_sens_s_latex = '$\ell_{I}$';
t_focus=1;
run('Figure_codes/Sensitivity_figures_in_time_s_Contour.m');



