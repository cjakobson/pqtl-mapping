function []=plot_hotspot_targets(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

    unique_genes=unique([all_pqtl_data.gene1; all_pqtl_data.gene2]);
    unique_genes=unique_genes(2:end);
    
    clear v1 v2
    for i=1:length(unique_genes)
        
        temp_idx=logical(ismember(all_pqtl_data.gene1,unique_genes{i})+...
            ismember(all_pqtl_data.gene2,unique_genes{i}));
        
        v1(i)=sum(temp_idx);
        %v2(i)=find(ismember(physicalGeneOrder,unique_genes{i}));
        
        %use median locus associated with the gene instead
        temp_idx=logical(ismember(variant_info.gene1,unique_genes{i})+...
            ismember(variant_info.gene2,unique_genes{i}));
        v2(i)=median(variant_info.index(temp_idx));
        
    end

    scatter(v2,v1,40,'k','filled',...
        'MarkerFaceAlpha',0.75)
    ylim([0 400])
    xlim([1 12054])
    set(gca,'YScale','log')

    

end