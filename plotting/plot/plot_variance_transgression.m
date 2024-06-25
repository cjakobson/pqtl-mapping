%plot heritability against mean abundance
function []=plot_variance_transgression(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    %calculate mapping stats to correlate
    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    for i=1:length(orf_names)

        temp_idx=ismember(all_pqtl_data.protein,orf_names{i});

        if sum(temp_idx>0)

            total_var_exp(i)=sum(all_pqtl_data.varExp(temp_idx));
            
        end

    end

    
    
    [transgressive_mat]=calculate_transgression(dependency_directory,output_directory);

    parent_max=max([transgressive_mat(:,1) transgressive_mat(:,2)],[],2);
    parent_min=min([transgressive_mat(:,1) transgressive_mat(:,2)],[],2);


    transgression=log2((transgressive_mat(:,4)./transgressive_mat(:,3)))-...
        log2((parent_max./parent_min));
    

    hold on
    axis square
    scatter(transgression,total_var_exp,10,'k','filled')
    xlim([-2 2])
    ylim([0 0.8])
    %set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    xlabel('transgression')
    ylabel('total var. exp.')
    
    
end

