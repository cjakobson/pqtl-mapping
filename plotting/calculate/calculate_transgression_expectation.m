function [p,var_ratio,var_labels,rm_mean,yjm_mean,f6_distr,f6_sim]=calculate_transgression_expectation(dependency_directory,output_directory)

    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    m=1;
    n=1;
    %change_thresh=0.1;
    change_thresh=Inf;
    p_thresh=10;
    rng(0)
    clear p
    %figure
    for i=1:length(orf_names)
    
        v_temp=input_mat(i,:);
    
        %RM and YJM means
        transgressive_mat(i,1)=mean(v_temp(rm_idx),'omitnan');
        transgressive_mat(i,2)=mean(v_temp(yjm_idx),'omitnan');
    
        parent_range(i)=abs(transgressive_mat(i,1)-transgressive_mat(i,2));
        %relative change
        parent_relative_change(i)=parent_range(i)/...
            mean([transgressive_mat(i,1) transgressive_mat(i,2)]);
    
        %compare to F6 for given relative change
        if parent_relative_change(i)<change_thresh
    
            f6_distr{m}=v_temp(f6_idx);
    
            %parent behavior
            temp_mean=mean([v_temp(rm_idx) v_temp(yjm_idx)],'omitnan');
            temp_std=std([v_temp(rm_idx) v_temp(yjm_idx)],[],'omitnan');
    
            %f6_distr{m}=(f6_distr{m}-temp_mean)/temp_std;
    
            n_sims=length(f6_idx);
    
            for j=1:n_sims
    
                temp_rand=randn;
    
                f6_sim{m}(j)=randn*temp_std+temp_mean;
    
            end
    
            %[h,p(m)]=kstest2(f6_distr{m},f6_sim{m});
            [h,p(m)]=vartest2(f6_distr{m},f6_sim{m});
            var_ratio(m)=var(f6_distr{m},[],'omitnan')/...
                var(f6_sim{m},[],'omitnan');
            var_labels{m}=orf_names{i};
            rm_mean(m)=mean(v_temp(rm_idx),'omitnan');
            yjm_mean(m)=mean(v_temp(yjm_idx),'omitnan');
    
            m=m+1;
    
        end
    
    end



end

