%plot replicates of RM against one another
function []=plot_reproducibility(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx]=...
        parse_raw_abundance(dependency_directory,output_directory);

    
    v1=log10(input_mat(:,rm_idx(1)));
    v2=log10(input_mat(:,rm_idx(2)));
    scatter(v1,v2,10,'filled','k','MarkerFaceAlpha',0.5)
    xlim([2 7])
    ylim(xlim)
    axis square
    [r p]=corr(v1,v2,'rows','complete');
    text(3,6,['r = ' num2str(r)])
    xlabel('RM11 rep1')
    ylabel('RM11 rep2')

    
end