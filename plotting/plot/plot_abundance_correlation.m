function [] = plot_abundance_correlation(gene_name1,gene_name2,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    protein_idx1=ismember(orf_names,gene_name1);
    protein_idx2=ismember(orf_names,gene_name2);
    
    v1=input_mat(protein_idx1,f6_idx);
    v2=input_mat(protein_idx2,f6_idx);

    hold on
    scatter(v1,v2,10,'k','filled')
    xlabel(gene_name1)
    ylabel(gene_name2)
    xlim([0 Inf])
    ylim([0 Inf])

    [r p]=corr(v1',v2','rows','complete');
    text(1e4,2e4,num2str(r))
    text(1e4,1e4,num2str(p))

    axis square


end


