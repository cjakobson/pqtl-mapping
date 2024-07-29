function []=plot_1K_mrna_protein(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    %include mapping prediction in boxplot
    
    
    cis_pqtn_data_mrna=calculate_1K_replication_mrna(dependency_directory,output_directory);
    
    cis_pqtn_data=calculate_1K_replication(dependency_directory,output_directory);
    
    v1=log2(cis_pqtn_data_mrna.rm_mean./cis_pqtn_data_mrna.yjm_mean);
    v2=log2(cis_pqtn_data.rm_mean./cis_pqtn_data.yjm_mean);
    
    hold on
    scatter(v1,v2,10,'k','filled')
    ylim([-1 1])
    xlim(ylim)
    axis square
    xlabel('mRNA effect (RM/YJM)')
    ylabel('protein effect (RM/YJM)')
    for i=1:height(cis_pqtn_data_mrna)
        
        if (abs(v1(i))>0.1)||((abs(v2(i))>0.1))
            
            text(v1(i),v2(i),cis_pqtn_data_mrna.commonName{i})
            
        end
        
    end
    
end


