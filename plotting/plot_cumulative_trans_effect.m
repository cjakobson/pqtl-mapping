function []= plot_cumulative_trans_effect(dependency_directory,output_directory)

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
    unique_proteins=unique(all_pqtl_data.protein);
    temp_idx1=all_pqtl_data.dist<1e3;
    for i=1:length(unique_proteins)
        
        temp_idx2=ismember(all_pqtl_data.protein,unique_proteins{i});
        
        %cis
        temp_idx=logical(temp_idx2.*temp_idx1);
        to_plot{1}(i)=sum(all_pqtl_data.varExp(temp_idx));
        
        %sum of trans
        temp_idx=logical(temp_idx2.*~temp_idx1);
        to_plot{2}(i)=sum(all_pqtl_data.varExp(temp_idx));
        
    end
    tempLabels={'sum of cis','sum of trans'};
    
    subplot(2,6,8)
    dot_and_line_plot(to_plot)
    ylim([0 0.15])
    xticks(1:length(to_plot))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel('mean varExp')
    
    
    mean(to_plot{1})
    mean(to_plot{2})
    
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(1.5,0.12,num2str(p))



end


