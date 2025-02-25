function [] = plot_pqtl_qtl_cumulative(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %pqtls
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    %growth qtls
    growth_data=readtable([dependency_directory 'linearNoRad.csv']);
    
    %cumulative varExp CDF
    clear to_plot
    unique_proteins=unique(all_pqtl_data.protein);
    for i=1:length(unique_proteins)
    
        temp_idx=find(ismember(all_pqtl_data.protein,unique_proteins{i}));
        [v_sorted,sort_idx]=sort(all_pqtl_data.varExp(temp_idx),'descend');
    
        for j=1:length(temp_idx)
    
            to_plot{1}(i,j)=sum(all_pqtl_data.varExp(temp_idx(sort_idx(1:j))));
    
        end
    
    end
    
    
    for i=1:height(growth_data)
        growth_data.merged_condition{i}=[growth_data.condition{i} '_' growth_data.time{i}];
    end
    
    unique_conditions=unique(growth_data.merged_condition);
    for i=1:length(unique_conditions)
        
        temp_idx=find(ismember(growth_data.merged_condition,unique_conditions{i}));
        [v_sorted,sort_idx]=sort(growth_data.varExp(temp_idx),'descend');
    
        for j=1:length(temp_idx)
    
            to_plot{2}(i,j)=sum(growth_data.varExp(temp_idx(sort_idx(1:j))));
    
        end
    
    end
    
    
    
    %subplot(2,4,3)
    hold on
    %pQTLs
    to_plot{1}(to_plot{1}==0)=nan;
    v1=median(to_plot{1},'omitnan');
    v2=std(to_plot{1},[],'omitnan');
    plot(v1,'k')
    plot(v1+v2,'--k')
    plot(v1-v2,'--k')
    %QTLs
    to_plot{2}(to_plot{2}==0)=nan;
    v1=median(to_plot{2},'omitnan');
    v2=std(to_plot{2},[],'omitnan');
    plot(v1,'b')
    plot(v1+v2,'--b')
    plot(v1-v2,'--b')
    ylim([0 0.6])
    xlim([0 25])
    axis square
    xlabel('(p)QTLs')
    ylabel('cumulative varExp')
    
    idx_to_test=10;
    [p h]=ranksum(to_plot{1}(:,idx_to_test),to_plot{2}(:,idx_to_test));
    text(idx_to_test,0.5,num2str(p))
    



end