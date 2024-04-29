function [] = plot_miami(gene_to_highlight,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %pQTLs
    pqtl_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    pqtl_input(pqtl_input.index==0,:)=[];   %geometric
    
    
    v1=pqtl_input.index;
    v2=pqtl_input.pVal;
    
    gene_idx1=logical(ismember(pqtl_input.gene1,gene_to_highlight)+...
        ismember(pqtl_input.gene2,gene_to_highlight));
    
    %growth QTLs
    growth_input=readtable([dependency_directory 'linearNoRad.csv']);
    growth_input(growth_input.index==0,:)=[];   %geometric
    
    
    v3=growth_input.index;
    v4=-growth_input.pVal;  %bottom side
    
    gene_idx2=logical(ismember(growth_input.gene1,gene_to_highlight)+...
        ismember(growth_input.gene2,gene_to_highlight));
    
    hold on
    scatter(v1,v2,v2,grey,'filled')
    scatter(v1(gene_idx1),v2(gene_idx1),v2(gene_idx1),'k','filled')
    
    scatter(v3,v4,-v4,grey,'filled')
    scatter(v3(gene_idx2),v4(gene_idx2),-v4(gene_idx2),'k','filled')
    
    ylim([-100 100])
    %just chromosome 15
    xlim([10010 11048])
    
end
    