function [transgressive_mat]=calculate_transgression(dependency_directory,output_directory)

    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
            
            
    %analyze transgressive inheritance
    for i=1:length(orf_names)

        v_temp=input_mat(i,:);

        %RM and YJM means
        transgressive_mat(i,1)=mean(v_temp(rm_idx),'omitnan');
        transgressive_mat(i,2)=mean(v_temp(yjm_idx),'omitnan');

        parent_range(i)=abs(transgressive_mat(i,1)-transgressive_mat(i,2));

        %f6 20th and 80th percentiles
        v_f6=v_temp(f6_idx);
        v_sort=sort(v_f6,'ascend');
        temp_length=length(v_sort);
        idx1=ceil(0.2*temp_length);
        idx2=ceil(0.8*temp_length);

        transgressive_mat(i,3)=v_sort(idx1);
        transgressive_mat(i,4)=v_sort(idx2);

        f6_range(i)=abs(transgressive_mat(i,3)-transgressive_mat(i,4));

    end
    
end


