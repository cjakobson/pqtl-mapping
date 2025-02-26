%plot heritability between parents
function []=plot_heritability_parents(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    [h_squared_rm,h_squared_yjm,h_squared_mean,~]=...
        calculate_heritability(dependency_directory,output_directory);
    

    hold on
    axis square
    scatter(h_squared_rm,h_squared_yjm,10,'k','filled')
    xlim([0 1])
    ylim(xlim)
    plot(xlim,ylim,':r')
    xlabel('H^2 RM11')
    ylabel('H^2 YJM975')
    
end

