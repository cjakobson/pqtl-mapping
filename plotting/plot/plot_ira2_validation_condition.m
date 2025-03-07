function [] = plot_ira2_validation_condition(condition_to_use,label_switch,dependency_directory,output_directory)

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

    [strain_volcano_p,strain_volcano_effect,mapping_beta,...
        rm_coherent,yjm_coherent,p_coherent,...
        condition_volcano_p,condition_volcano_effect,n_diff_exp]=calculate_condition_stats(dependency_directory,output_directory);

    condition_idx=find(ismember(condition_names,condition_to_use));

    v1=strain_volcano_effect{condition_idx};
    v2=mapping_beta;
    hold on
    scatter(log2(v1),v2,10,'k','filled')
    axis square
    xlim([-1 1])
    ylim([-0.6 0.6])
    xlabel('log_2 YJM975 IRA2*/WT')
    ylabel('mapping \beta')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')
    title(condition_names(condition_idx))
    
    
    
    if label_switch
        for i=1:length(protein_names)
            text(log2(v1(i)),v2(i),protein_names(i))
        end
    end




end