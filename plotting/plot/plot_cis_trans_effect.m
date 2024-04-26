function []= plot_cis_trans_effect(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];

    clear to_plot
    %all global pQTLs
    to_plot{1}=all_pqtl_data.varExp;


    temp_idx1=~isinf(all_pqtl_data.dist);
    %diff chr
    to_plot{2}=all_pqtl_data.varExp(~temp_idx1);
    %same chr
    to_plot{3}=all_pqtl_data.varExp(temp_idx1);
    
    temp_idx2=all_pqtl_data.dist<1e3;
    %outside 1kb, same chr
    to_plot{4}=all_pqtl_data.varExp(logical(temp_idx1.*~temp_idx2));
    %within 1kb
    to_plot{5}=all_pqtl_data.varExp(temp_idx2);
    
    temp_labels={'all','diff. chr','same chr','same chr >1kB','same chr <=1kB'};
    
    dot_and_line_plot(to_plot)
    ylim([0 0.06])
    xticks(1:length(to_plot))
    xtickangle(45)
    xticklabels(temp_labels)
    ylabel('mean varExp')
    
    [p h]=ranksum(to_plot{1},to_plot{5});
    text(3,0.05,num2str(p))
    
    
    


end


