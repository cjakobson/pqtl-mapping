 function [] = plot_qq(protein_idx,dependency_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    [true_p_sorted, expected_p_sorted] = ...
        calculate_qq(protein_idx,dependency_directory);
    
    scatter(expected_p_sorted,true_p_sorted)
    
    
end

