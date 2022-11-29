% This script plots how a variable sens_var_s moves in time t, together with the succession of projections in time s
figure
x= (1:T+S+1);
index=1;
for sens_run_2=1:sens_run_2
    for sens_run_1=1:sens_run_1
    f=getfield(sens_res(sens_run_1,sens_run_2),sens_var_s);
    %c(index,:)=f(1,:);
    plot(x,f)
    color = rand(1,3);
    hold on
    xlabel('time t')
    ylabel(sens_var_s)
    index=index+1;
    end
end