function []= plot_structure_all_possible(property_idx,dependency_directory,output_directory)

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
    %ts only
    temp_idx1=mis_ts==1;
    to_plot{2}=properties_sim{1,property_idx}(temp_idx1);
    %tv only
    to_plot{3}=properties_sim{1,property_idx}(~temp_idx1);
    %all segregating
    to_plot{4}=properties_all_segregating{property_idx};
    
    hold on
    easy_box(to_plot)
    ylim([0 yLim2(property_idx)])
    xticklabels({'all poss.','Ts','Tv','all segr.'})
    text(4,10,num2str(median(to_plot{4},'omitnan')/median(to_plot{1},'omitnan')))
    % for i=1:3
    %     [p h]=ranksum(to_plot{i},to_plot{4});
    %     text((i+4)/2,35,num2str(p))
    % end
    title(property_labels{property_idx})


end


