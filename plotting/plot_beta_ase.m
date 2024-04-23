%plot beta vs ASE
function []=plot_beta_ase(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    all_pqtl_data(all_pqtl_data.dist>1e3,:)=[];

    qtn_idx=all_pqtl_data.isQtn==1;

    hold on
    %v1=all_pqtl_data.majorAfNoRad;
    v1=all_pqtl_data.aseTagRmAfNoRad;
    v1(isnan(v1))=all_pqtl_data.maxRmAfNoRad(isnan(v1));
    %v1(all_pqtl_data.maxRmAfNoRad<all_pqtl_data.majorAfNoRad)=1-v1(all_pqtl_data.maxRmAfNoRad<all_pqtl_data.majorAfNoRad);
    v2=all_pqtl_data.beta;
    scatter(v1,v2,10,'k','filled')
    scatter(v1(qtn_idx),v2(qtn_idx),20,'r','filled')
    xlabel('RM11 AF no rad')
    ylabel('\beta')
    xlim([0 1])
    ylim([-1 1])
    hold on
    plot(xlim,[0 0],':r')
    plot([0.5 0.5],ylim,':r')
    axis square

    
end


