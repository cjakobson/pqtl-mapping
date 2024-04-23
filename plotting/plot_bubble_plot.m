%plot beta vs ASE
function []=plot_bubble_plot(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    p_ratio=3;
    beta_ratio=75;
    
    all_pqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    all_pqtl(all_pqtl.bPos==1,:)=[];
    
    [physical_gene_order,breaks_to_plot,breaks_to_plot2]=...
        calculate_chr_breaks(dependency_directory,output_directory);

    %put dot plot in here too
    hold on
    xlim([1 12054])
    ylim([1 length(physical_gene_order)])
    %axis square
    for i=1:length(breaks_to_plot)
        plot([breaks_to_plot(i) breaks_to_plot(i)],ylim,'k','LineWidth',1)
    end
    for i=1:length(breaks_to_plot2)
        plot(xlim,[breaks_to_plot2(i) breaks_to_plot2(i)],'k','LineWidth',1)
    end
    for i=1:length(physical_gene_order)
        %all pQTLs
        temp_idx=ismember(all_pqtl.protein,physical_gene_order{i});
        temp_pqtl=all_pqtl.index(temp_idx);
        temp_p=all_pqtl.pVal(temp_idx);
        temp_beta=all_pqtl.beta(temp_idx);
        rm_idx=temp_beta>0;
        yjm_idx=temp_beta<0;
        scatter(temp_pqtl(rm_idx),i*ones(1,length(temp_pqtl(rm_idx))),p_ratio*abs(temp_p(rm_idx)),...
            'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.5,...
            'MarkerEdgeColor',blue,'LineWidth',1)
        scatter(temp_pqtl(yjm_idx),i*ones(1,length(temp_pqtl(yjm_idx))),p_ratio*abs(temp_p(yjm_idx)),...
            'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.5,...
            'MarkerEdgeColor',orange,'LineWidth',1)
    end
    xlabel('genomic coordinate')
    ylabel('ORFs in genome order')
    xlim([1 12054])
    ylim([1 length(physical_gene_order)])
    axis off
    
end


