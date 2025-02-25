%rarefaction plots by descending abundance
function []=plot_marginal_rarefaction(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
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

    %random permutations of protein order
    for j=1:100
    
        rng(j)
    
        sort_idx=randperm(length(orf_names));
    
        for i=1:length(sort_idx)
    
            temp_query=orf_names(sort_idx(1:i));
        
            v_temp=all_pqtl.bestCandidate(logical(ismember(all_pqtl.protein,temp_query).*...
                (all_pqtl.dist>1e3)));
        
            mat_to_plot(j,i)=length(unique(v_temp));
    
        end
    
    end
    
    v_random=median(mat_to_plot);
    
    %plot 20/80%
    for i=1:length(v_random)
        
        [~,temp_sort_idx]=sort(mat_to_plot(:,i));
    
        temp_idx1=floor(0.2*length(temp_sort_idx));
        temp_idx2=ceil(0.8*length(temp_sort_idx));
    
        v_20(i)=mat_to_plot(temp_sort_idx(temp_idx1),i);
        v_80(i)=mat_to_plot(temp_sort_idx(temp_idx2),i);
    
    end



    %plot marginal pQTLs
    v1=movmean(v_to_plot(2:end)-v_to_plot(1:(end-1)),25);
    v2=movmean(v_random(2:end)-v_random(1:(end-1)),25);
    
    
    %subplot(2,4,4)
    hold on
    plot(1:length(v1),v1,'k')
    %plot(1:length(v_to_plot2),v_to_plot2)
    plot(1:length(v2),v2,'b')
    axis square
    xlim([0 length(v1)])
    xlabel('number of proteins measured')
    ylabel('number of marginal pQTLs')
    title('rarefaction')
    ylim([0 10])
    
    temp_split=floor(length(v2)/5);
    mean1=mean(v2(1:temp_split));
    mean2=mean(v2((end-temp_split):end));
    text(500,4,['1st vs 5th quintile yields ' num2str(mean1/mean2) '-fold'])


    
end



