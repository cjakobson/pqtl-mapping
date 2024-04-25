function [v_mean,v_background,p_val,sign_sum1,sign_sum2] = calculate_sign_test(dependency_directory,output_directory)


    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];
    

    %slide pVal thresh and ask about coherence
    v_bins=1:20;
    
    v_protein=unique(all_pqtl_data.protein);
    
    f_mat=nan(length(v_protein),length(v_bins));
    max_mat=f_mat;
    min_mat=f_mat;
    sign_mat=nan(length(v_protein),length(v_bins));
    clear to_plot v_mean v_background
    for i=1:length(v_bins)
    
        temp_table=all_pqtl_data(all_pqtl_data.pVal>v_bins(i),:);
    
        temp_signs=temp_table.beta;
        
        temp_sum1=sum(temp_signs>0);
        temp_sum2=sum(temp_signs<0);
    
        temp_max=max([temp_sum1 temp_sum2]);
        temp_min=min([temp_sum1 temp_sum2]);
    
        v_background(i)=temp_max/(temp_max+temp_min);
    
        for j=1:length(v_protein)
    
            temp_idx=ismember(temp_table.protein,v_protein{j});
    
            if sum(temp_idx)>1
        
                temp_signs=temp_table.beta(temp_idx);
        
                temp_sum1=sum(temp_signs>0);
                temp_sum2=sum(temp_signs<0);
        
                temp_max=max([temp_sum1 temp_sum2]);
                temp_min=min([temp_sum1 temp_sum2]);
        
                f_mat(j,i)=temp_max/(temp_max+temp_min);
                
                max_mat(j,i)=temp_max;
                min_mat(j,i)=temp_min;
    
                if temp_sum1>temp_sum2
                    sign_mat(j,i)=1;
                elseif temp_sum1<temp_sum2
                    sign_mat(j,i)=-1;
                end
    
            end
    
        end
    
        to_plot{i}=f_mat(:,i);
        v_mean(i)=mean(to_plot{i},'omitnan');
        v_N(i)=sum(~isnan(to_plot{i}));
        sign_sum1(i)=sum(sign_mat(:,i)==1,'omitnan');
        sign_sum2(i)=sum(sign_mat(:,i)==-1,'omitnan');
        
    end
    
    v_max=sum(max_mat,'omitnan');
    v_min=sum(min_mat,'omitnan');
    
    
    %pVal vs background
    for i=1:length(v_max)
        
        p_val(i)=binocdf(v_max(i)/(v_min(i)+v_max(i)),v_min(i)+v_max(i),v_background(i));
        
    end





end