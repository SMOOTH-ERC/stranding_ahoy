% This script plots the movement of var_s in time s, evaluated at a specific period t_focus in time t
hold off
x= (1:T+S+1);
index=1;
for sens_run_2=1:sens_run_2
    for sens_run_1=1:sens_run_1
        f=getfield(sens_res(sens_run_1,sens_run_2),var_s);
        c(index,:)=f(t_focus,:);
        color = rand(1,3);
        plot(x,c(index,:))
        hold on
        index=index+1;
    end
end
legend(var_s,'Fontsize',16)
xlabel("Future Time Periods From Present")
ylabel("Expected Utilisation Rate")
title('Expected Utilisation Rates for Differing Expected Transition Speeds');
