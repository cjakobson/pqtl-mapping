%plot heritability against mean abundance
function []=plot_heritability_abundance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    
    [~,~,h_squared_mean,~]=...
        calculate_heritability(dependency_directory,output_directory);
    

    hold on
    axis square
    scatter(h_squared_mean,mean(input_mat,2,'omitnan'),10,'k','filled')
    xlim([0 1])
    ylim([1e3 2e6])
    set(gca,'YScale','log')
    %plot(xlim,ylim,':r')
    xlabel('mean H^2')
    ylabel('mean abundance')
    
    
end

