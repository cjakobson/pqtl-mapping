function []= plot_structure_rare_common(property_idx,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    af_thresh=1e-3;


    [properties_1K,properties_sim,properties_all_segregating,properties_all_other,properties_pqtn,...
        struct_mis_1K,mis_secondary,mis_af] = parse_structure_analysis(dependency_directory,output_directory);

    property_labels={'accessible surface area','neighbors'};

    
    structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
    mis_idx=mis_af<af_thresh;
    
    %subset by structure
    clear to_plot
    clear temp_labels
    
    m=1;
    for j=1:length(structure_labels)
        
        %observed
        temp_idx1=struct_mis_1K'==j;
        temp_idx2=mis_idx;
        
        to_plot(m)=mean(properties_1K{1,property_idx}(logical(temp_idx1.*~temp_idx2)),'omitnan')/...
            mean(properties_1K{1,property_idx}(logical(temp_idx1.*temp_idx2)),'omitnan');
        temp_labels{m}=structure_labels{j};
        [p_val(m) h]=ranksum(properties_1K{1,property_idx}(logical(temp_idx1.*~temp_idx2)),...
            properties_1K{1,property_idx}(logical(temp_idx1.*temp_idx2)));
        m=m+1;
        
    end
    
    hold on
    bar(to_plot,'BaseValue',1)
    xticks(1:length(temp_labels))
    xtickangle(45)
    xticklabels(temp_labels)
    ylabel([property_labels{property_idx} ' common/rare'])
    xlim([0.5 length(temp_labels)+0.5])
    %ylim([2/3 3/2])
    ylim([0.85 1.2])
    set(gca,'YScale','log')
    title('missense 1K')
    %p_vals
    for j=1:length(p_val)
        text(j,1.15,num2str(p_val(j)))
    end
    
    
    
    
end


