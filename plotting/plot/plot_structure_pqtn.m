function []= plot_structure_pqtn(property_idx,dependency_directory,output_directory)

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
    to_plot{2}=properties_all_other{property_idx};
    to_plot{3}=properties_pqtn{property_idx};
    
    hold on
    easy_box(to_plot)
    ylim([0 yLim2(property_idx)])
    xticklabels({'all poss.','all other','pQTNs'})
    title(property_labels{property_idx})

    [p h]=ranksum(to_plot{2},to_plot{3});
    text(2.5,35,num2str(p))
    
    
end


