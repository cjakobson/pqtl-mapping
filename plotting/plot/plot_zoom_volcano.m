function []=plot_zoom_volcano(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %to get protein names
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);

    [fold_change,p_val]=calculate_parental_mean_fc(dependency_directory,output_directory);
    
    pqtl_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    %don't count OD600
    pqtl_input=pqtl_input(pqtl_input.index>0,:);

    for i=1:length(orf_names)
        
        n_pqtls(i)=sum(ismember(pqtl_input.protein,orf_names{i}));

    end

    sum(n_pqtls>0)

    sum(p_val>0.05)

    sum((n_pqtls>0).*(p_val'>0.05))

    mean(n_pqtls(p_val'>0.05))
    
    has_pqtl_idx=n_pqtls>0;
    
    hold on   
    scatter(log2(fold_change(has_pqtl_idx)),-log10(p_val(has_pqtl_idx)),25,'k','filled')
    scatter(log2(fold_change(~has_pqtl_idx)),-log10(p_val(~has_pqtl_idx)),25,'or')

    ylim([0 -log10(0.05)])
    xlim([-0.4 0.4])
    
    
    %xticks(-2:5)

%     plot_idx=p_val<0.05;
%     rm_idx=fold_change>1;
%     yjm_idx=fold_change<1;
%     
%     scatter(log2(fold_change(logical(plot_idx.*rm_idx))),-log10(p_val(logical(plot_idx.*rm_idx))),25,blue,'filled')%,'MarkerFaceAlpha',0.5)
%     scatter(log2(fold_change(logical(plot_idx.*yjm_idx))),-log10(p_val(logical(plot_idx.*yjm_idx))),25,orange,'filled')%,'MarkerFaceAlpha',0.5)

    ylabel('-log_{10}q')
    xlabel('log_2FC')
    axis square

end


    