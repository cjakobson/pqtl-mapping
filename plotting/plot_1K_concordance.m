function []=plot_1K_concordance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    cis_pqtn_data=calculate_1K_replication(dependency_directory,output_directory);
    
    

    hold on
    v1=cis_pqtn_data.beta;
    v2=log2(cis_pqtn_data.rm_mean./cis_pqtn_data.yjm_mean);
    to_plot{1}=v2(v1<0);
    to_plot{2}=v2(v1>0);
    easy_box(to_plot)
    ylim([-0.4 0.4])
    xticks(1:2)
    xticklabels({'\beta<0','\beta>0'})
    ylabel('log_2RM allele mean/YJM allele mean')
    title('all local pQTLs')
    %[h p]=ttest2(to_plot{1},to_plot{2});
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(1.5,-0.3,['p = ' num2str(p)])



    
    
end