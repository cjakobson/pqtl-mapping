function [] = plot_pqtn_qtn_scores_overlap(plot_offset,gene_name,locus_to_use,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    

    
    [pqtn_sorted, qtn_sorted] = ...
        calculate_pqtn_qtn_scores(gene_name,locus_to_use,dependency_directory,output_directory);


    subplot(2,4,plot_offset+1)
    hold on
    for i=1:length(pqtn_sorted)
        plot((1:length(pqtn_sorted{i})),pqtn_sorted{i},'k','LineWidth',1)
        title(num2str(locus_to_use))
        [temp_max,max_idx]=max(pqtn_sorted{i});
        scatter(max_idx,temp_max,25,'k','filled')

    end
    %xlim([-9 Inf])
    %plot([11 11],ylim,':r')
    axis square


%     subplot(2,4,plot_offset+2)
%     hold on
    for i=1:length(qtn_sorted)
        plot((1:length(qtn_sorted{i})),qtn_sorted{i},'k','LineWidth',0.5)
        %title(num2str(locus_to_use(j)-10))
        [temp_max,max_idx]=max(qtn_sorted{i});
        scatter(max_idx,temp_max,25,'k','filled')

    end
    xlim([0 19])
    %plot([11 11],ylim,':r')
    axis square
    
    
    
end

