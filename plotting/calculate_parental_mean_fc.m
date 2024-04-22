function [fold_change,p_val]=calculate_parental_mean_fc(dependency_directory,output_directory)

    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names]=...
        parse_raw_abundance(dependency_directory,output_directory);

    rm_mean=mean(input_mat(:,rm_idx),2,'omitnan');
    yjm_mean=mean(input_mat(:,yjm_idx),2,'omitnan');
    


    %FC
    fold_change=rm_mean./yjm_mean;

    %p_val
    rawp_val=zeros(size(fold_change));
    for i=1:length(orf_names)

        v1=input_mat(i,rm_idx);
        v2=input_mat(i,yjm_idx);

        [~, rawp_val(i)]=ttest2(v1,v2);

    end
    %bonferroni
    %p_val=p_val.*(length(orf_names));
    %b-h
    p_val=mafdr(rawp_val,'BHFDR',true);


end

