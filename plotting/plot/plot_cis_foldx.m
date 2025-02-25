function [] = plot_cis_foldx(dependency_directory,output_directory)

    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    %missense variants - destabilizing => lower abundance?
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    
    %cis missense pQTNs
    all_pqtl_data=all_pqtl_data(all_pqtl_data.isQtn==1,:);
    all_pqtl_data=all_pqtl_data(ismember(all_pqtl_data.variantType,'missense_variant'),:);
    all_pqtl_data=all_pqtl_data(all_pqtl_data.dist<1e3,:);
    
    %set rm == ref
    v_beta=all_pqtl_data.beta;
    v_beta(all_pqtl_data.refAllele==all_pqtl_data.yjmAllele)=...
        -v_beta(all_pqtl_data.refAllele==all_pqtl_data.yjmAllele);
    
    v_foldx=all_pqtl_data.variantFoldX;
    
    
    foldx_thresh=1;
    
    clear to_plot
    to_plot{1}=v_beta(v_foldx<foldx_thresh);
    to_plot{2}=v_beta(v_foldx>foldx_thresh);
    
    
    %subplot(2,4,8)
    hold on
    scatter(v_foldx,-v_beta,10,'k','filled')
    axis square
    ylim([-0.6 0.6])
    xlim([-5 15])
    ylabel('pQTL effect (ref-alt)')
    xlabel('FoldX')
    plot([1 1],ylim,':r')
    plot(xlim,[0 0],':r')
    
    



end


