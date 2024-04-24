function []=plot_transgression(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
        
    
    [transgressive_mat]=calculate_transgression(dependency_directory,output_directory);

    %sort by grand mean
    [~,sort_idx]=sort(mean(input_mat,2,'omitnan'),'ascend');

    hold on
    scatter(1:length(orf_names),transgressive_mat(sort_idx,1),25,'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.5)
    scatter(1:length(orf_names),transgressive_mat(sort_idx,2),25,'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.5)

    scatter(1:length(orf_names),transgressive_mat(sort_idx,3),25,'filled','k','MarkerFaceAlpha',0.5)
    scatter(1:length(orf_names),transgressive_mat(sort_idx,4),25,'filled','k','MarkerFaceAlpha',0.5)

    set(gca, 'YScale', 'log')
    xlim([0 length(orf_names)])
    xlabel('proteins (sorted by mean abundance)')
    ylabel('abundance (AU)')
    legend({'RM11','YJM975','F_6 20th percentile','F_6 80th percentile'},...
        'Location','Northwest')


        
        
end

