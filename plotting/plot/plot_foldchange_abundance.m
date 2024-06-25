%plot heritability against mean abundance
function []=plot_foldchange_abundance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    
    [fold_change,p_val]=calculate_parental_mean_fc(dependency_directory,output_directory);

    fold_change=log2(fold_change);
    

    hold on
    axis square
    scatter(mean(input_mat,2,'omitnan'),fold_change,10,'k','filled')
    ylim([-2 5])
    xlim([1e3 2e6])
    set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    ylabel('log_2fold change')
    xlabel('mean abundance')
    
    
end

