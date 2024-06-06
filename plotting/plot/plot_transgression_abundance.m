%plot heritability against mean abundance
function []=plot_transgression_abundance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    
    [transgressive_mat]=calculate_transgression(dependency_directory,output_directory);

    parent_max=max([transgressive_mat(:,1) transgressive_mat(:,2)],[],2);
    parent_min=min([transgressive_mat(:,1) transgressive_mat(:,2)],[],2);


    transgression=log2((transgressive_mat(:,4)./transgressive_mat(:,3)))-...
            log2((parent_max./parent_min));
    

    hold on
    axis square
    scatter(mean(input_mat,2,'omitnan'),transgression,10,'k','filled')
    ylim([-3 3])
    xlim([1e3 2e6])
    set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    ylabel('mean H^2')
    xlabel('mean abundance')
    
    
end

