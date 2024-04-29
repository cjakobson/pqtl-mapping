function []= plot_structure_histogram(property_idx,increment,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    af_thresh=1e-3;


    [properties_1K,properties_sim,properties_all_segregating,properties_all_other,properties_pqtn,...
        struct_mis_1K,mis_secondary,mis_af,mis_ts] = parse_structure_analysis(dependency_directory,output_directory);

    property_labels={'accessible surface area','neighbors'};

    tic
    
    yLim1=[0 0];
    yLim2=[250 40];

    %compare all segregating to all possible
    clear to_plot
    
    %missense
    to_plot{1}=properties_sim{1,property_idx};
    to_plot{2}=properties_all_segregating{property_idx};
    
    hold on
    histogram(to_plot{1},0:increment:yLim2(property_idx),'Normalization','Probability')
    histogram(to_plot{2},0:increment:yLim2(property_idx),'Normalization','Probability')
    xlim([0 yLim2(property_idx)])
    ylim([0 Inf])
    title(property_labels{property_idx})
    axis square
    legend({'sim.','all segr.'})

end


