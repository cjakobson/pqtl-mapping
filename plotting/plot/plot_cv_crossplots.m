
function []=plot_cv_crossplots(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    rm_mat=input_mat(:,rm_idx);
    yjm_mat=input_mat(:,yjm_idx);
    f6_mat=input_mat(:,f6_idx);
    sgrp_mat=input_mat(:,sgrp_idx);


    for i=1:length(orf_names)

        rm_cv(i)=std(rm_mat(i,:),[],'omitnan')/mean(rm_mat(i,:),'omitnan');
        yjm_cv(i)=std(yjm_mat(i,:),[],'omitnan')/mean(yjm_mat(i,:),'omitnan');

        f6_cv(i)=std(f6_mat(i,:),[],'omitnan')/mean(f6_mat(i,:),'omitnan');
        sgrp_cv(i)=std(sgrp_mat(i,:),[],'omitnan')/mean(sgrp_mat(i,:),'omitnan');
        
        %parent_cv(i)=mean([temp_rm_cv temp_yjm_cv]);

    end

    to_plot{1}=rm_cv;
    to_plot{2}=yjm_cv;
    to_plot{3}=f6_cv;

    axis_labels={'RM CV','YJM CV','F_6 CV'};
   
    m=1;
    for i=1:length(to_plot)

        for j=(i+1):length(to_plot)

            v1=to_plot{i};
            v2=to_plot{j};
            subplot(2,4,m)
            hold on
            axis square
            scatter(v2,v1,10,'k','filled')
            ylim([0 1])
            xlim([0 1])
            %set(gca,'XScale','log')
            %plot(xlim,ylim,':r')
            xlabel(axis_labels{i})
            ylabel(axis_labels{j})

    
            [r p]=corr(v1',v2','rows','complete');
            v_xlim=xlim;
            v_ylim=ylim;
            text(0.8*v_xlim(2),0.9*v_ylim(2),['r = ' num2str(r)])
            text(0.8*v_xlim(2),0.8*v_ylim(2),['p = ' num2str(p)])
            

            m=m+1;

        end

    end


    %plot parents and segregants norm to grand median; sort by
    %transgression
    [p,var_ratio,var_labels,rm_mean,yjm_mean,f6_distr,f6_sim]=...
        calculate_transgression_expectation(dependency_directory,output_directory);

    [v_sort,sort_idx]=sort(var_ratio);

    for i=1:length(var_ratio)
    
        %temp_median=median(f6_distr{sort_idx(i)},'omitnan');

        v1(i)=rm_mean(sort_idx(i));%/temp_median;
        v2(i)=yjm_mean(sort_idx(i));%/temp_median;
    
        % v_temp=sort(f6_distr{sort_idx(i)},'ascend');
        % v3(i)=v_temp(floor(0.2*length(v_temp)))/temp_median;
        % v4(i)=v_temp(ceil(0.8*length(v_temp)))/temp_median;

        temp_mean=mean(f6_distr{sort_idx(i)},'omitnan');
        v3(i)=abs(temp_mean-v1(i))/v1(i);
        v4(i)=abs(temp_mean-v2(i))/v2(i);

    end
    
    subplot(2,4,4)
    hold on
    % for i=1:length(v3)
    %     plot([v3(i) v4(i)],[i i],'k')
    % end
    % scatter(v1,1:length(v1),10,blue,'filled')
    % scatter(v2,1:length(v1),10,orange,'filled')
    scatter(v3,v4,10,'k','filled')
    axis square
    % set(gca,'xscale','log')
    % set(gca,'yscale','log')
    xlim([0 1])
    ylim(xlim)
    xlabel('|(F_6 mean-RM mean)/RM mean|')
    ylabel('|(F_6 mean-YJM mean)/YJM mean|')
    plot(xlim,ylim,':r')
    text(0.2,0.9,num2str(sum(v4>v3)))
    text(0.9,0.2,num2str(sum(v4<v3)))
    
end

