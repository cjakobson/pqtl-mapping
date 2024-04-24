function []=plot_npqtls_per_protein(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);


    %nQtls per protein
    clear to_plot
    unique_proteins=unique(all_pqtl_data.protein);
    for i=1:length(unique_proteins)
        
        temp_idx=ismember(all_pqtl_data.protein,unique_proteins{i});
        
        to_plot(i)=sum(temp_idx);
        
    end
    
    
    histogram(to_plot)
    xlim([0 40])
    ylim([0 150])
    xlabel('number of pQTLs')
    ylabel('number of proteins')
    axis square


end