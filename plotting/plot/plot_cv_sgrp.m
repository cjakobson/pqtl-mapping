%plot heritability against mean abundance
function []=plot_cv_sgrp(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    rm_mat=input_mat(:,rm_idx);
    yjm_mat=input_mat(:,yjm_idx);
    f6_mat=input_mat(:,f6_idx);
    sgrp_mat=input_mat(:,sgrp_idx);


    for i=1:length(orf_names)

        temp_rm_cv=std(rm_mat(i,:),[],'omitnan')/mean(rm_mat(i,:),'omitnan');
        temp_yjm_cv=std(yjm_mat(i,:),[],'omitnan')/mean(yjm_mat(i,:),'omitnan');

        f6_cv(i)=std(f6_mat(i,:),[],'omitnan')/mean(f6_mat(i,:),'omitnan');
        sgrp_cv(i)=std(sgrp_mat(i,:),[],'omitnan')/mean(sgrp_mat(i,:),'omitnan');
        
        parent_cv(i)=mean([temp_rm_cv temp_yjm_cv]);

    end

   
    v1=f6_cv./parent_cv;
    v2=sgrp_cv./parent_cv;
    

    hold on
    axis square
    scatter(v2,v1,10,'k','filled')
    ylim([0 5])
    xlim([0 10])
    %set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    xlabel('C.V. SGRP/C.V. parents')
    ylabel('C.V. F6/C.V. parents')
    
    
    [r p]=corr(v1',v2','rows','complete');
    v_xlim=xlim;
    v_ylim=ylim;
    text(0.8*v_xlim(2),0.9*v_ylim(2),num2str(r))
    text(0.8*v_xlim(2),0.8*v_ylim(2),num2str(p))
    
    
end

