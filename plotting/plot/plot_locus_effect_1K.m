function []=plot_locus_effect_1K(gene_name,common_name,locus_number,y_lim,...
    dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    cis_pqtn_data=calculate_1K_replication(dependency_directory,output_directory);
    
    temp_idx=find(ismember(cis_pqtn_data.protein,gene_name));
    to_plot{1}=cis_pqtn_data.rm_array{temp_idx};
    to_plot{2}=cis_pqtn_data.yjm_array{temp_idx};
    
    hold on
    easy_box(to_plot)
    title(cis_pqtn_data.commonName{temp_idx})
    xticks(1:2)
    xtickangle(45)
    xticklabels({'RM allele','YJM allele'})
    ylim([0 y_lim])
    v_temp=ylim;
    ylim([0 v_temp(2)])
    text(1,0.8*v_temp(2),['\beta = ' num2str(cis_pqtn_data.beta(temp_idx))])
    for j=1:2
        text(j,0.1*v_temp(2),num2str(length(to_plot{j})))
    end
    
    
end


