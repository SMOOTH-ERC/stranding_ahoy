% Displays trend over time in 3D with 1 variable sens_par_2 changing, with
% sens.res.??? as the z axis. 
% ell_I_res = [sens_res.ell_I];
% [X,Y] = meshgrid(1:T+1, sens_par_range_2);
% mesh(Y, X, transpose(ell_I_res));
% ylabel("Time");
% xlabel("Planning Horizon");
% zlabel("Share of Investment in Renewables");
% title('Trend Of Renewable Investment for Differing Planning Horizons')

% Figure plotting var_t (to be chosen below as a script) in time t for each loop
% Should only have one sensitivity variable
% var_t='ell_I';
% run('Figure_codes\Figures_in_time_t');

% Figure plotting var_s (to be chosen below as a script) in time s, evaluated at a 
% specific period t_focus (to be chosen below, usually =1) in time t
%  var_s='sigma_u_H';
%  t_focus=1;
%  run('Figure_codes/Figures_in_time_s.m');

% Figure plotting var_sens_s (to be chosen below as a script) against sensitivity
% parameters, evaluated at a specific period t_focus (to be chosen below, usually =1) in time t
% Returns nothing is no sensitivity range is specified
% Returns a 2D plot if there is one sensitivity parameter
% Returns a 3D plot if there are two sensitivity parameters
% var_sens_s1 = 'ell_I';
% var_sens_s2 = 'change_K';
% var_sens_s3 = 'change_E';
% t_focus1=2;
% t_focus2=2;
% run('Figure_codes/Sensitivity_figures_in_time_s_multiple');
% 
%var_sens_s = 'ell_I_alt';
%var_sens_s_latex = strcat('$\alpha=$', sprintf('%.1f',alpha), "");
%t_focus=1;
%run('Figure_codes/Sensitivity_figures_in_time_s');

%var_sens_s = 'ell_I_alt2';
%var_sens_s_latex = '$\ell_{I,alt2}$';
%t_focus=1;
%run('Figure_codes/Sensitivity_figures_in_time_s');

% Figure plotting how sens_var_s (to be chosen below as a script) moves in
% both time t and time s
% sens_var_s  =   'mu_pi_H';
% run('Figure_codes\Figures_in_time_t_and_s');
% 
%  plot([1:T+1], transpose(u_H), "r", [1:T+1], min(transpose(mu_u_H)), "b", [1:T+S+1], transpose(mu_u_H(1,:)), "g--", [1:T+S+1], transpose(mu_u_H(20,:)), "g--",   [1:T+S+1], transpose(mu_u_H(40,:)), "g--");%, transpose(mu_u_H(30,:)), "g--", [1:71], transpose(mu_u_H(40,:)), "g--");
%  set(gca, 'FontSize', 10);
%  xlabel('Time');
%  ylabel('Utilisation Rate');
%  title('Utilisation and Stranding');
%  legend("Actual utilisation rate", "Maximum Mean Expected Stranding", "Expected Stranding (t = 1)", "Expected Stranding (t = 20)", "Expected Stranding (t = 40)");%, "Expected Stranding (t = 30)", "Expected Stranding (t = 40)", 'Location', "SouthEast", 'Fontsize', 8);

%plot([1:T+1], transpose(u_H), "r", [1:T+1], min(transpose(mu_u_H)), "b");
%  
% plot([1:T+1], transpose(ell_I), "r", [1:T+1], transpose(ell_K), "g", [1:T+1], transpose(ell_E), "b");
% set(gca, 'FontSize', 10);
% xlabel('Time');
% ylabel('Share');
% title('Key Energy Ratios');
% legend("Investment", "Capital", "Energy", 'Location', "SouthEast");

%plot([1:71], transpose(mu_u_H(10,:)), "r");

%plot([1:S], ell_E_e(1,1:S));
%set(gca, 'FontSize', 20);
%xlabel('Time', 'interpreter', 'latex');
%ylabel('Share', 'interpreter', 'latex');
%title('Expected share of low-carbon energy in total demand $E(\ell_{E})$', 'interpreter', 'latex');


%% Cleaning
% Delete several parameters to improve readability of results
%clear('beta','eps_A', 'g_e', 'hb_uH', 'lb_uH','psi_H', 'psi_L', 'rho', 's', 'S', 'sens_par_1', 'sens_par_2', 'sens_par_name_1', 'sens_par_name_2', 'sens_par_range_1', 'sens_par_range_2', 'sens_run_1', 'sens_run_2', 'sens_var_s', 'T', 't_focus', 't_s', 'var_s', 'var_sens_s')


% Some functions for refence here: not used at the moment
function out = getVarName(var)
    out = inputname(1);
end  

function out = sens(struct)
    out = struct(1,1).var
end


%pd1 = makedist('Stable','alpha',1,'beta',-1,'gam',0.02,'delta',0.95);
%pd2 = makedist('Stable','alpha',1,'beta',-1,'gam',0.02,'delta',0.95);
%pd3 = makedist('Stable','alpha',1,'beta',-1,'gam',0.02,'delta',0.95);



%pdf2 = pdf(pd2,x);
%pdf3 = pdf(pd3,x);


%x = -1:0.01:1.5;
%pdf1 = pdf(pd1,x);
%figure
%plot(x,pdf1,'b-');
%hold on
%plot(x,pdf2,'r-.');
%plot(x,pdf3,'k--');
%title('Compare Alpha Parameters in Stable Distribution PDF Plots')
%legend('\alpha = 2','\alpha = 1','\alpha = 0.5','Location','northwest')
%hold off

%plot(1:20, sens_res(sens_run_1,sens_run_2).ell_E_e);
