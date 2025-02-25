function []= plot_sign_test_jackknife(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    n_splits=10;
    
    hold on
    for i=1:n_splits
        rng(i)
        [v_mean,v_background,p_val,sign_sum1,sign_sum2] = ...
                calculate_sign_test_split(dependency_directory,output_directory);
            
        plot(v_mean,'-k')
        plot(v_background,'--k')
    end
    ylim([0.4 1])
    plot(xlim,[0.5 0.5],':r')
    xlabel('p threshold')
    ylabel('f_{coherent}')
    axis square



    

end


