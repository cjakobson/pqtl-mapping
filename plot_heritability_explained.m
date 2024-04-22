%plot heritability against mean abundance
function []=plot_heritability_explained(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [~,~,h_squared_mean,var_exp]=...
        calculate_heritability(dependency_directory,output_directory);
    

    hold on
    to_plot=var_exp./h_squared_mean;
    to_plot(to_plot==0)=nan;
    histogram(to_plot,0:0.05:1)
    ylabel('frequency')
    xlabel('fraction of H^2 explained')
    xlim([0 1])
    axis square
    
    
end

