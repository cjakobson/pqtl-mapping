function []=plot_pqtn_blosum(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];
    
    
    variant_info=readtable([dependency_directory  'variantInfoStructure.csv']);
    
    
    qtn_idx=all_pqtl_data.isQtn==1;
    missense_idx=ismember(all_pqtl_data.variantType,'missense_variant');
    trans_idx=all_pqtl_data.dist>1e3;
    to_plot{1}=all_pqtl_data.variantFoldX(logical(qtn_idx.*missense_idx.*trans_idx));
    %exclude missing
    to_plot{1}(to_plot{1}==0)=[];
    
    other_idx=~ismember(variant_info.index,all_pqtl_data.index(logical(qtn_idx.*missense_idx.*trans_idx)));
    missense_idx=ismember(variant_info.variantType,'missense_variant');
    to_plot{2}=variant_info.variantFoldX(logical(missense_idx.*other_idx));
    to_plot{2}(to_plot{2}==0)=[];
    
    hold on
    easy_box(to_plot)
    ylim([-1 5])
    xticklabels({'pQTNs','all other'})
    xtickangle(45)
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(1.5,3,num2str(p))
    



end


