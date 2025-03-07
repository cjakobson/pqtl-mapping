function [] = plot_condition_diff_exp(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [mean_mat,v_label_mean,v_strain_mean,v_condition_mean,...
        clean_mat,v_label_all,v_strain_all,v_condition_all,...
        z_mat,protein_names]=parse_ira2_conditions(dependency_directory,output_directory);

    condition_names=unique(v_condition_all);
    strain_names=unique(v_strain_all);

    mapping_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    [strain_volcano_p,strain_volcano_effect,mapping_beta,...
        rm_coherent,yjm_coherent,p_coherent,...
        condition_volcano_p,condition_volcano_effect,n_diff_exp]=calculate_condition_stats(dependency_directory,output_directory);

    short_names={'EtOH','flc','glc','mal','rap','teb'};
    
    [v_sort,sort_idx]=sort(n_diff_exp,'ascend');


    hold on
    bar(n_diff_exp(sort_idx)./length(protein_names))
    xticks(1:length(short_names))
    xticklabels(short_names(sort_idx))
    axis square
    ylim([0 1])
    xlim([0 length(short_names)+1])
    title('proteome remodeling')
    ylabel('fraction diff. exp.')
    xlabel('condition')





    

end