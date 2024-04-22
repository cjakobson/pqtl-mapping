function []=plot_sensitivity_simulations(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [f_discovered,v_bins,temp_labels]=...
        calculate_sensitivity_simulations(dependency_directory,output_directory);
    
    hold on
    plot(1-f_discovered)
    ylim([0 1])
    xlim([0 length(v_bins)])
    xticks(1:length(temp_labels))
    xticklabels(temp_labels)
    xtickangle(45)
    axis square
    ylabel('fraction discovered')
    xlabel('effect size')
    
    

end


