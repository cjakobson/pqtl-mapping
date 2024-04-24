%plot protein levels in parents and F6
function []=plot_parents_f6_boxplot(gene_to_plot,dependency_directory,output_directory)

    %gene_to_plot is systematic name of ORF

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    input_data=readtable([dependency_directory '211031_SegregantProteomicsData_DetectionThreshold80_genes_ORF.tsv'],...
        'FileType','text');
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    rm_cols=ismember(strain_names,'RM11');
    yjm_cols=ismember(strain_names,'YJM975');

    f6_cols=ismember(strain_names,'F6');
    
    protein_idx=ismember(input_data.Protein_Group,gene_to_plot);
    
    to_plot{1}=input_mat(protein_idx,rm_cols);
    to_plot{2}=input_mat(protein_idx,yjm_cols);
    to_plot{3}=input_mat(protein_idx,f6_cols);
    
    easy_box_with_dots(to_plot)
    title(gene_to_plot)
    ylim([0 Inf])
    
    
    
end