
function []=plot_validation2(gene_name,mutant_to_plot,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    %input_data=readtable([dependency_directory 'Segregants_validation2_annotated.xlsx']);
    input_data=readtable([dependency_directory '20240607_proteomics-data_validation_set2.xlsx']);
    
    strain_names=unique(input_data.strain_id);
    protein_names=unique(input_data.COMMON);
    protein_names(cellfun(@isempty,protein_names))=[];

    
    protein_idx=logical(ismember(input_data.COMMON,gene_name));

    
    strain_idx1=ismember(input_data.strain_id,'YDJ6635');
    tempIdx=logical(protein_idx.*strain_idx1);
    to_plot{1}=input_data.abundance(tempIdx);
    
    strain_idx2=ismember(input_data.strain_id,'YDJ6649');
    tempIdx=logical(protein_idx.*strain_idx2);
    to_plot{2}=input_data.abundance(tempIdx);

    strain_idx3=ismember(input_data.strain_id,mutant_to_plot);
    tempIdx=logical(protein_idx.*strain_idx3);
    to_plot{3}=input_data.abundance(tempIdx);
    
    
    hold on
    for j=1:length(to_plot)
        v_mean(j)=mean(to_plot{j});
    end
    bar(v_mean)
    for j=1:length(to_plot)
        %length(to_plot{j})
        scatter(j*ones(length(to_plot{j}),1),to_plot{j},25,'k','filled')
    end
    xticks(1:3)
    xtickangle(45)
    xticklabels({'YJM975','RM11',mutant_to_plot})
    title(gene_name)
    [h p]=ttest2(to_plot{2},to_plot{3});
    plot([2 3],[1.1*mean(v_mean(2:3)) 1.1*mean(v_mean(2:3))],'-k')
    text(2.5,1.2*mean(v_mean(2:3)),num2str(p))

    
end


