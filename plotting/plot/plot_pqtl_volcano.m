function []=plot_pqtl_volcano(gene_name,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    pqtl_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    trans_idx=pqtl_input.dist>1e3;
    od_idx=pqtl_input.index==0;

    pqtl_to_use=pqtl_input(logical(trans_idx.*~od_idx),:);


    %quantify specifically for FRE1
    pqtl_idx=find(ismember(pqtl_to_use.gene1,gene_name)+...
        ismember(pqtl_to_use.gene2,gene_name));


    %pseudovolcano from mapping
    v1=pqtl_to_use.beta(pqtl_idx);
    v2=pqtl_to_use.pVal(pqtl_idx);

    hold on
    scatter(v1(v1>0),v2(v1>0),50,blue,'filled')
    scatter(v1(v1<0),v2(v1<0),50,orange,'filled')
    axis square
    xlabel('\beta')
    ylabel('pVal')
    ylim([0 40])
    xlim([-0.3 0.3])
    
end
