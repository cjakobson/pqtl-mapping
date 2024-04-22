function []=plot_od_boxplot(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    %final OD boxplot
    load([dependency_directory 'rmOd.mat'])
    load([dependency_directory 'yjmOd.mat'])


    to_plot{1}=rmOd;
    to_plot{2}=yjmOd;

    easy_box_with_dots(to_plot)
    xticklabels({'RM','YJM'})
    ylim([0 2])
    
end

