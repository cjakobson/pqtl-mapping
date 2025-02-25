function [] = plot_pqtl_qtl_complexity(dependency_directory,output_directory)


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
    

    %number of loci
    clear to_plot
    unique_proteins=unique(all_pqtl_data.protein);
    for i=1:length(unique_proteins)
        to_plot{1}(i)=sum(ismember(all_pqtl_data.protein,unique_proteins{i}));
    end
    
    
    for i=1:height(growth_data)
        growth_data.merged_condition{i}=[growth_data.condition{i} '_' growth_data.time{i}];
    end
    
    unique_conditions=unique(growth_data.merged_condition);
    for i=1:length(unique_conditions)
        to_plot{2}(i)=sum(ismember(growth_data.merged_condition,unique_conditions{i}));
    end
    
    %subplot(2,8,2)
    hold on
    easy_box(to_plot)
    xticklabels({'pQTLs','QTLs'})
    set(gca,'YScale','log')
    ylabel('number of (p)QTLs')
    ylim([1 1e3])
    xlim([0.5 2.5])
    for i=1:length(to_plot)
        text(i,1,num2str(median(to_plot{i})),'Rotation',45)
    end
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(1.5,700,num2str(p))
    



end