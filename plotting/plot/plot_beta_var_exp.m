
function []=plot_beta_var_exp(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    all_pqtl_data(all_pqtl_data.index==0,:)=[];   %geometric
    
    
    v1=abs(all_pqtl_data.beta);
    v2=all_pqtl_data.varExp;
    scatter(v1,v2,5,'k','filled')
    axis square
    set(gca,'YScale','log')
    xlim([0 0.5])
    xlabel('\beta')
    ylabel('variance explained')

    
end



