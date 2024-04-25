function [] = plot_pqtl_maf(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

    mapping_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv'],...
        'FileType','text');

    v1=variant_info.nMinor;
    v2=mapping_data.nMinor(mapping_data.index>0);
    
    hold on
    histogram(v1,0:20:500,'Normalization','Probability')
    histogram(v2,0:20:500,'Normalization','Probability')
    axis square
    xlabel('minor allele frequency')
    ylabel('norm. freq.')

end