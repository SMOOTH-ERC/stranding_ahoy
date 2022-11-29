% This script plots how a variable var_sens_s changes for a range of parameter values (sens_par_range_1 and
% sens_par_range_2), evaluated at a specific period t_focus in time t

% If no sensitivity ranges are specified, nothing happens
if  length(sens_par_range_1)==1 && length(sens_par_range_2)==1
    % do nothing
% If one range of values is specified, return a 2D sensitivity plot for the parameter whose range is specified
elseif length(sens_par_range_1)==1 || length(sens_par_range_2)==1
    if length(sens_par_range_1)==1 
        sens_par_range  =   sens_par_range_2;
        sens_par_name   =   sens_par_name_2;
        sens_run        =   sens_run_2;
    else
        sens_par_range  =   sens_par_range_1;
        sens_par_name   =   sens_par_name_1;
        sens_run        =   sens_run_1;
    end
    index=1;
  fig =  figure;
    for sens_run=1:sens_run
        f=getfield(sens_res(sens_run),var_sens_s);
        var_f(index,:)=f(t_focus,:);
        color = rand(1,3);
        index=index+1;
    end
    plot(sens_par_range,var_f)
    ylabel(var_sens_s)
    xlabel(sens_par_name)
% If two ranges of values are specified, return a 3D plot
else
   fig = figure('WindowState','maximized', 'Position',[0.130520833333333,0.172827225130891,0.775000000000001,0.815000000000001]);
    [x_axis,y_axis]     = meshgrid(sens_par_range_1,sens_par_range_2);
    for  sens_run_2=1:sens_run_2
        for sens_run_1=1:sens_run_1
            f=getfield(sens_res(sens_run_1,sens_run_2),var_sens_s);
            z_axis(sens_run_2,sens_run_1) = f(t_focus,:); 
            z_axis(isnan(z_axis)) = 0;
        end
    end
        legend_name =var_sens_s;
        [C,h] = contour(x_axis,y_axis,z_axis, "ShowText","off", LineWidth=3);
        %clabel(C,h,'FontSize',15);
        contourLegend(h);
        ax.Fontsize = 50;
        xlabel(sens_par_name_1,'Fontsize',45, 'interpreter', 'latex', 'Rotation',0);
        ylabel(sens_par_name_2,'Fontsize',45, 'interpreter', 'latex', 'Rotation',90);
        zlabel(sens_var_name,'Fontsize',45, 'interpreter', 'latex');
        xl = get(gca, 'XLabel');
        Fontsize = 45;
        xlFontsize = get(xl, 'Fontsize');
        xAx = get(gca, 'XAxis');
        set(xAx, 'Fontsize', Fontsize);
        set(xl, 'Fontsize', xlFontsize);
        yl = get(gca, 'YLabel');
        ylFontsize = get(yl, 'Fontsize');
        yAx = get(gca, 'YAxis');
        set(yAx, 'Fontsize', Fontsize);
        set(yl, 'Fontsize', ylFontsize);
        zl = get(gca, 'ZLabel');
        zlFontsize = get(zl, 'Fontsize');
        zAx = get(gca, 'ZAxis');
        set(zAx, 'Fontsize', Fontsize);
        set(zl, 'Fontsize', zlFontsize);
        %title(['Value of ', var_sens_s_latex ,' at time 1 depending on parameter values'], 'interpreter',  'latex', 'Fontsize', 17)
        %title(var_sens_s_latex, 'interpreter',  'latex', 'Fontsize', 35)
     %ylabel('ell_I');
     ax.Position = [0.130520833333333,0.172827225130891,0.775000000000001,0.815000000000001];
end 