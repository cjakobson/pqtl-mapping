function [] = plot_pqtl_microarray(gene_name,deletion_name,dependency_directory,output_directory)

    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    array_input=readtable([dependency_directory 'deleteome_all_mutants_controls.txt']);

    idx_to_use=find(ismember(all_pqtl_data.gene1,gene_name)+...
        ismember(all_pqtl_data.gene2,gene_name));
    
    clear mapping_beta protein_label array_z_score
    for i=1:length(idx_to_use)
    
        temp_locus=gene_name;%all_pqtl_data.gene1{idxToUse(i)};
        temp_target=all_pqtl_data.protein{idx_to_use(i)};
    
        mapping_beta(i)=all_pqtl_data.beta(idx_to_use(i));
        protein_label{i}=all_pqtl_data.protein{idx_to_use(i)};
    
        %temp_idx1=ismember(deletionGeneSystematic,tempLocus);
        temp_idx1=find(ismember(array_input.Properties.VariableNames,...
            [deletion_name '_delVs_Wt']));
        
        temp_idx2=find(ismember(array_input.systematicName,temp_target));
    
        if (~isempty(temp_idx1))&&(~isempty(temp_idx2))
    
            v_temp=table2array(array_input(temp_idx2,temp_idx1));
            array_z_score(i)=v_temp;
            
        else
            
            array_z_score(i)=nan;
    
        end
    
    end
    
    
    hold on
    scatter(array_z_score,mapping_beta,10,'k','filled')
    axis square
    xlim([-3 3])
    ylim([-0.8 0.8])
    xlabel('microarray vs WT')
    ylabel('mapping \beta')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')
    title(gene_name)
    
    [r p]=corr(array_z_score',mapping_beta','rows','complete');
    text(2,0.4,num2str(r))
    text(2,0.3,num2str(p))
    

end