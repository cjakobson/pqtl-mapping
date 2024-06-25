%plot heritability against mean abundance
function []=plot_pqtls_foldchange(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    %calculate mapping stats to correlate
    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    
    for i=1:length(orf_names)

        temp_idx=ismember(all_pqtl_data.protein,orf_names{i});

        n_pqtls(i)=sum(temp_idx);

    end

    
    rm_mat=input_mat(:,rm_idx);
    yjm_mat=input_mat(:,yjm_idx);
    f6_mat=input_mat(:,f6_idx);

    for i=1:length(orf_names)

        temp_rm_cv=std(rm_mat(i,:),[],'omitnan')/mean(rm_mat(i,:),'omitnan');
        temp_yjm_cv=std(yjm_mat(i,:),[],'omitnan')/mean(yjm_mat(i,:),'omitnan');

        cv1(i)=mean([temp_rm_cv temp_yjm_cv]);

        cv2(i)=std(f6_mat(i,:),[],'omitnan')/mean(f6_mat(i,:),'omitnan');

    end

    cv1=cv1';
    cv2=cv2';

    norm_cv=cv2./cv1;

    

    hold on
    axis square
    scatter(norm_cv,n_pqtls,10,'k','filled')
    xlim([0 4])
    ylim([0 35])
    %set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    xlabel('norm. CV')
    ylabel('n pQTLs')
    
    [r p]=corr(norm_cv,n_pqtls','rows','complete');
    v_xlim=xlim;
    v_ylim=ylim;
    text(0.8*v_xlim(2),0.9*v_ylim(2),num2str(r))
    text(0.8*v_xlim(2),0.8*v_ylim(2),num2str(p))
    
    
end

