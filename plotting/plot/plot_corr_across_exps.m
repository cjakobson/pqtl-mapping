function []=plot_corr_across_exps(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [prot_corr, prot_corr_1K, prot_corr_5K] =...
        calculate_correlations(dependency_directory,output_directory);
    
    v1=reshape(prot_corr,[],1);
    v2=reshape(prot_corr_1K,[],1);
    
    subplot(2,4,1)
    hold on
    hist3([v1 v2],'Ctrs',{-0.6:0.02:0.6, -0.6:0.02:0.6},'CdataMode','auto')%,'EdgeColor','none','FaceColor','interp')
    xlim([-0.6 0.6])
    ylim(xlim)
    axis square
    xlabel('F6')
    ylabel('1K')
    view(2)
    
    
    
    
    v1=reshape(prot_corr,[],1);
    v2=reshape(prot_corr_5K,[],1);
    
    subplot(2,4,2)
    hold on
    hist3([v1 v2],'Ctrs',{-0.6:0.02:0.6, -0.6:0.02:0.6},'CdataMode','auto')%,'EdgeColor','none','FaceColor','interp')
    xlim([-0.6 0.6])
    ylim(xlim)
    axis square
    xlabel('F6')
    ylabel('5K')
    view(2)
    
    
    

end

