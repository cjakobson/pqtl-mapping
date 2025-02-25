%rarefaction plots by descending abundance
function []=plot_transgression_example(gene_to_plot,dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [p,var_ratio,var_labels,rm_mean,yjm_mean,f6_distr,f6_sim]=...
        calculate_transgression_expectation(dependency_directory,output_directory);
    
    idx_to_use=find(ismember(var_labels,gene_to_plot));
        
    hold on
    f6_max=max([f6_distr{idx_to_use} f6_sim{idx_to_use}],[],'omitnan');
    n_bins=20;
    histogram(f6_distr{idx_to_use},0:f6_max/n_bins:f6_max,'Normalization','probability')
    histogram(f6_sim{idx_to_use},0:f6_max/n_bins:f6_max,'Normalization','probability')
    axis square
    title(var_labels{idx_to_use})
    xlabel('SWATH-MS abundance (A.U.)')
    ylabel('rel. freq.')
    text(0,0.1,num2str([p(idx_to_use)]))
    ylim([0 0.35])
    plot([rm_mean(idx_to_use) rm_mean(idx_to_use)],ylim,'Color',blue)
    plot([yjm_mean(idx_to_use) yjm_mean(idx_to_use)],ylim,'Color',orange)
    legend({'actual F_6','sim F_6','RM','YJM'})


    
    
end



