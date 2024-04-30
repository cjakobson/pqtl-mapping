function [] = plot_pqtl_qtl_effect_size(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %pqtls
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    all_pqtl_data(all_pqtl_data.index==0,:)=[];   %geometric
    
    %growth qtls
    growth_data=readtable([dependency_directory 'linearNoRad.csv']);
    growth_data(growth_data.index==0,:)=[];   %geometric
    

    to_plot{1}=all_pqtl_data.varExp;
    to_plot{2}=growth_data.varExp;
    
    hold on
    easy_box(to_plot)
    xticklabels({'pQTLs','QTLs'})
    set(gca,'YScale','log')
    ylim([1e-5 1])
    ylabel('var. exp.')


end