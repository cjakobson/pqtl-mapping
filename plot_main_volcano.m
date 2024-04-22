function []=plot_main_volcano(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [fold_change,p_val]=calculate_parental_mean_fc(dependency_directory,output_directory);
    
    
    hold on
    scatter(log2(fold_change),-log10(p_val),5,'k','filled')

    ylim([0 60])
    xlim([-2 5])
    xticks(-2:5)

    plot_idx=p_val<0.05;
    rm_idx=fold_change>1;
    yjm_idx=fold_change<1;
    
    scatter(log2(fold_change(logical(plot_idx.*rm_idx))),-log10(p_val(logical(plot_idx.*rm_idx))),25,blue,'filled')%,'MarkerFaceAlpha',0.5)
    scatter(log2(fold_change(logical(plot_idx.*yjm_idx))),-log10(p_val(logical(plot_idx.*yjm_idx))),25,orange,'filled')%,'MarkerFaceAlpha',0.5)

    ylabel('-log_{10}q')
    xlabel('log_2FC')
    axis square

end


    