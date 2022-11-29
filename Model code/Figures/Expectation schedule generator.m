figure
plot(1:30, sens_res(1,sens_run_2).ell_E_e(1,1:30));
hold on
for j = 2:7
plot(1:30, sens_res(j,sens_run_2).ell_E_e(1,1:30));
end
%legend('$b_{\ell} = 0.1$','$b_{\ell} = 0.15$','$b_{\ell} = 0.2$','$b_{\ell} = 0.25$', '$b_{\ell} = 0.3$', '$b_{\ell} = 0.35$', '$b_{\ell} = 0.4$', 'interpreter', 'latex', 'location', "northwest")
xlabel('Time $s$','interpreter', 'latex')
ylabel('$E(\ell_{E}$)','interpreter', 'latex')
hold off

figure
plot(1:30, sens_res(1,sens_run_2).mu_u_H(1,1:30));
hold on
for j = 2:7
plot(1:30, sens_res(j,sens_run_2).mu_u_H(1,1:30));
end
legend('$b_{\ell} = 0.1$','$b_{\ell} = 0.15$','$b_{\ell} = 0.2$','$b_{\ell} = 0.25$', '$b_{\ell} = 0.3$', '$b_{\ell} = 0.35$', '$b_{\ell} = 0.4$', 'interpreter', 'latex', 'location', "southwest")
xlabel('Time $s$','interpreter', 'latex')
ylabel('$E(u^{*}_{H}$)','interpreter', 'latex')
hold off



