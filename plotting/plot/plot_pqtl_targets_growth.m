function [] = plot_pqtl_targets_growth(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    

    %how often are hotspots QTNs? are they enriched?
    
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    all_pqtl_data(all_pqtl_data.index==0,:)=[];   %geometric
    %all_pqtl_data=all_pqtl_data(all_pqtl_data.isQtn==1,:);
    
    pqtl_genes=unique([all_pqtl_data.gene1; all_pqtl_data.gene2]);
    pqtl_genes=pqtl_genes(2:end);
    
    
    growth_data=readtable([dependency_directory 'linearNoRad.csv']);
    growth_data(growth_data.index==0,:)=[];   %geometric
    growth_data=growth_data(growth_data.isQtn==1,:);


    clear n_targets
    for i=1:length(pqtl_genes)
        
        temp_idx=logical(ismember(all_pqtl_data.gene1,pqtl_genes{i})+...
            ismember(all_pqtl_data.gene2,pqtl_genes{i}));
        
        n_targets(i)=sum(temp_idx);
    
        temp_idx=logical(ismember(growth_data.gene1,pqtl_genes{i})+...
            ismember(growth_data.gene2,pqtl_genes{i}));
    
        n_growth(i)=length(unique(growth_data.condition(temp_idx)));
        is_growth(i)=logical(sum(temp_idx));
        
    end
    
    
    %include all genes w SNPs as a control
    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);
    
    all_genes=unique([variant_info.gene1; variant_info.gene2]);
    clear is_growth_all
    for i=1:length(all_genes)
    
        temp_idx=logical(ismember(growth_data.gene1,all_genes{i})+...
            ismember(growth_data.gene2,all_genes{i}));
    
        is_growth_all(i)=logical(sum(temp_idx));
    
    end
    
    
    
    clear v1 v2
    target_thresh=1:15;
    for i=1:length(target_thresh)
    
        temp_idx=n_targets>=target_thresh(i);
    
        v1(i)=mean(is_growth(temp_idx));
        v2(i)=std(is_growth(temp_idx))/sqrt(sum(temp_idx));
    
    end
    
    
    
    
    
    %subplot(2,4,6)
    hold on
    errorbar(target_thresh,v1,v2,'.k')
    axis square
    ylabel('fraction impacting growth')
    xlabel('no. of targets')
    xlim([0 max(target_thresh)])
    ylim([0.25 0.6])
    plot(xlim,[mean(is_growth) mean(is_growth)],':k')
    plot(xlim,[mean(is_growth_all) mean(is_growth_all)],':r')
    title('pQTL genes')
    



end