function []=plot_od_correlations(gene_name,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    load([dependency_directory 'pQTLfilename.mat'])
    load([dependency_directory 'pQTLtrait.mat'])
    load([dependency_directory 'f6od.mat'])
    
    temp_idx=find(ismember(filename,gene_name));
    
    v1=f6od;
    v2=trait{temp_idx};
    
    scatter(v1,v2,10,'k','filled')
    axis square
    title(gene_name)
    xlim([0 2.5])
    ylim([-4 4])
    xlabel('OD')
    ylabel('normalized abundance')

end


    