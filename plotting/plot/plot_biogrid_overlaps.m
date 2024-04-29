function []=plot_biogrid_overlaps(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [bar_mat] = calculate_biogrid_overlaps(dependency_directory,output_directory);
    
    bar(bar_mat,'stacked')
    ylabel('fraction')
    title('BioGrid')
    
end