%plot heritability against mean abundance
function []=plot_heritability_cv(dependency_directory,output_directory)

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

    for i=1:length(orf_names)

        temp_rm_cv=std(rm_mat(i,:),[],'omitnan')/mean(rm_mat(i,:),'omitnan');
        temp_yjm_cv=std(yjm_mat(i,:),[],'omitnan')/mean(yjm_mat(i,:),'omitnan');

        cv(i)=mean([temp_rm_cv temp_yjm_cv]);

    end

    cv=cv';

    
    [~,~,h_squared_mean,~]=...
        calculate_heritability(dependency_directory,output_directory);
    

    hold on
    axis square
    scatter(cv,h_squared_mean,10,'k','filled')
    ylim([0.2 1])
    xlim([0 0.65])
    %set(gca,'XScale','log')
    %plot(xlim,ylim,':r')
    ylabel('mean H^2')
    xlabel('mean C.V.')
    
    
end

