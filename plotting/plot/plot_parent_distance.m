function []=plot_parent_distance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    

    %plot parents and segregants norm to grand median; sort by
    %transgression
    [p,var_ratio,var_labels,rm_mean,yjm_mean,f6_distr,f6_sim]=...
        calculate_transgression_expectation(dependency_directory,output_directory);

    [v_sort,sort_idx]=sort(var_ratio);

    for i=1:length(var_ratio)

        v1(i)=rm_mean(sort_idx(i));%/temp_median;
        v2(i)=yjm_mean(sort_idx(i));%/temp_median;

        temp_mean=mean(f6_distr{sort_idx(i)},'omitnan');
        v3(i)=abs(temp_mean-v1(i))/v1(i);
        v4(i)=abs(temp_mean-v2(i))/v2(i);

    end
    

    hold on
    scatter(log2(v3),log2(v4),10,'k','filled')
    axis square
    xlim([-8 0])
    ylim(xlim)
    xlabel('|(F_6 mean-RM mean)/RM mean|')
    ylabel('|(F_6 mean-YJM mean)/YJM mean|')
    plot(xlim,ylim,':r')
    text(-7,0,num2str(sum(v4>v3)))
    text(0,-7,num2str(sum(v4<v3)))
    
end

