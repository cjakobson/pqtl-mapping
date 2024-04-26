function []= plot_sign_test_by_strain(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;


    [v_mean,v_background,p_val,sign_sum1,sign_sum2] = ...
        calculate_sign_test(dependency_directory,output_directory);
    
    hold on
    plot(sign_sum1)
    plot(sign_sum2)
    legend({'RM','YJM'})
    xlabel('p threshold')
    ylabel('N_{coherent proteins}')
    axis square


    

end


