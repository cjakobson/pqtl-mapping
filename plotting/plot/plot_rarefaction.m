%rarefaction plots by descending abundance
function []=plot_rarefaction(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    
    
    v_mean=mean(input_mat,2,'omitnan');
    v_mean(isnan(v_mean))=-Inf;
    
    
    
    all_pqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    downsample_pqtl=readtable([dependency_directory 'linearPqtlOdDownsample_FDR_0.1.csv']);
    
    

    [v_sorted,sort_idx]=sort(v_mean,'descend');

    clear v_to_plot v_to_plot2 v_to_plot3
    for i=1:length(sort_idx)

        temp_query=orf_names(sort_idx(1:i));

        v_temp=all_pqtl.bestCandidate(logical(ismember(all_pqtl.protein,temp_query).*...
            (all_pqtl.dist>1e3)));

        v_to_plot(i)=length(unique(v_temp));

        v_temp=all_pqtl.bestCandidate(logical(ismember(all_pqtl.protein,temp_query).*...
            (all_pqtl.dist<=1e3)));

        v_to_plot2(i)=length(unique(v_temp));


        v_temp=downsample_pqtl.bestCandidate(logical(ismember(downsample_pqtl.protein,temp_query).*...
            (downsample_pqtl.dist>1e3)));

        v_to_plot3(i)=length(unique(v_temp));

    end


    
    hold on
    plot(1:length(v_to_plot),v_to_plot)
    plot(1:length(v_to_plot2),v_to_plot2)
    plot(1:length(v_to_plot3),v_to_plot3)
    axis square
    xlabel('number of proteins measured')
    ylabel('number of unique pQTLs')
    title('by descending abunance')
    xlim([0 length(v_to_plot)])
    ylim([0 2500])
    
end



