function []= plot_structure_ts_tv(property_idx,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    af_thresh=1e-3;

    
    yLim1=[0 0];
    yLim2=[250 40];

    [properties_1K,properties_sim,properties_all_segregating,properties_all_other,properties_pqtn,...
        struct_mis_1K,mis_secondary,mis_af,mis_ts] = parse_structure_analysis(dependency_directory,output_directory);

    property_labels={'accessible surface area','neighbors'};

    
    structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
    clear to_plot temp_labels pVal
    m=1;
    for j=1:length(structure_labels)
        
        %observed
        temp_idx1=struct_mis_1K'==j;
        
        %sim
        temp_idx2=mis_secondary==j;
        
        %all observed
        to_plot{m}=properties_1K{1,property_idx}(temp_idx1);
        temp_labels{m}=[structure_labels{j} '1K'];
        m=m+1;
        
        
        %all sim
        to_plot{m}=properties_sim{1,property_idx}(temp_idx2);
        temp_labels{m}=[structure_labels{j} 'all sim'];
        m=m+1;
        
        
        temp_idx3=mis_ts==1;
        %sim Ts
        to_plot{m}=properties_sim{1,property_idx}(logical(temp_idx2.*temp_idx3));
        temp_labels{m}=[structure_labels{j} 'Ts sim'];
        m=m+1;
        
        %sim Tv
        to_plot{m}=properties_sim{1,property_idx}(logical(temp_idx2.*~temp_idx3));
        temp_labels{m}=[structure_labels{j} 'Tv sim'];
        m=m+1;
        
        
    end
    
    hold on
    %bar(to_plot,'BaseValue',1)
    easy_box(to_plot)
    xticks(1:length(temp_labels))
    xtickangle(45)
    xticklabels(temp_labels)
    ylabel(property_labels{property_idx})
    title('missense')
    %xlim([0.5 length(temp_labels)+0.5])
    %ylim([0.8 1.3])
    ylim([0 yLim2(property_idx)])
    %set(gca,'YScale','log')
    
    
    
end


