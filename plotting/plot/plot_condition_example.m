function [] = plot_condition_example(protein_to_use,condition_to_use,y_lim,...
    dependency_directory,output_directory)

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

    temp_protein_idx=ismember(protein_names,protein_to_use);
    conditions_to_plot={'Glucose',condition_to_use};

    m=1;
    for j=1:length(conditions_to_plot)
    
        temp_condition_idx=ismember(v_condition_all,conditions_to_plot{j});

        temp_strain_idx=ismember(v_strain_all,'YJM');
        to_plot{m}=clean_mat(temp_protein_idx,...
            logical(temp_condition_idx.*temp_strain_idx));
        m=m+1;

        temp_strain_idx=ismember(v_strain_all,'YJM-IRA2-AG');
        to_plot{m}=clean_mat(temp_protein_idx,...
            logical(temp_condition_idx.*temp_strain_idx));
        m=m+1;

    end

    %subplot(2,8,4+i)
    easy_box_with_dots(to_plot)
    ylabel('abundance (AU)')
    xticks(1:4)
    xticklabels({'WT','edit','WT','edit'})
    title(protein_to_use)
    ylim([0 y_lim])
    xlim([0 5])


end