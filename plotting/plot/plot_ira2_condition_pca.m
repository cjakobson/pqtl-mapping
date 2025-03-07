function [] = plot_ira2_condition_pca(dependency_directory,output_directory)

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

    [coeff,score,latent,tsquared,explained,mu]=pca(z_mat');

    %plot PCA
    v_condition_color={'k','r','g','b','c','m','y'};
    v_strain_symbol={'o','o'};
    
    
    hold on
    
    v1=score(:,1);
    v2=score(:,2);
    m=1;
    for i=1:length(condition_names)
    
        condition_idx=ismember(v_condition_all,condition_names{i});
    
        for j=1:length(strain_names)
    
            strain_idx=ismember(v_strain_all,strain_names{j});
    
            temp_idx=logical(condition_idx.*strain_idx);
    
            if j==1
                scatter(v1(temp_idx),v2(temp_idx),25,v_strain_symbol{j},...
                    'MarkerFaceColor',v_condition_color{i},...
                    'MarkerEdgeColor','none')
            elseif j==2
                scatter(v1(temp_idx),v2(temp_idx),25,v_strain_symbol{j},...
                    'MarkerEdgeColor',v_condition_color{i},...
                    'MarkerFaceColor','none')

            end
    
            legend_entries{m}=[condition_names{i} '-' strain_names{j}];
            m=m+1;
    
        end
    
    end
    
    
    axis square
    legend(legend_entries)
    
end