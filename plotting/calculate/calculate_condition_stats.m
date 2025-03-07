function [strain_volcano_p,strain_volcano_effect,mapping_beta,...
    rm_coherent,yjm_coherent,p_coherent,...
    condition_volcano_p,condition_volcano_effect,n_diff_exp]=calculate_condition_stats(dependency_directory,output_directory)

    [mean_mat,v_label_mean,v_strain_mean,v_condition_mean,...
        clean_mat,v_label_all,v_strain_all,v_condition_all,...
        z_mat,protein_names]=parse_ira2_conditions(dependency_directory,output_directory);

    condition_names=unique(v_condition_all);
    strain_names=unique(v_strain_all);

    mapping_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);



    for i=1:length(condition_names)
    
        wt_idx=ismember(v_strain_all,'YJM');
        edit_idx=ismember(v_strain_all,'YJM-IRA2-AG');
    
        condition_idx=ismember(v_condition_all,condition_names{i});
    
        strain_volcano_p{i}=nan(length(protein_names),1);
        strain_volcano_effect{i}=nan(length(protein_names),1);  %also use to compare to mapping
        
        for j=1:length(protein_names)
    
            temp_idx1=logical(condition_idx.*edit_idx);
            v_temp1=clean_mat(j,temp_idx1);
    
            temp_idx2=logical(condition_idx.*wt_idx);
            v_temp2=clean_mat(j,temp_idx2);
    
            [h,strain_volcano_p{i}(j)]=ttest2(v_temp1,v_temp2);
            strain_volcano_effect{i}(j)=mean(v_temp1,'omitnan')/...
                mean(v_temp2,'omitnan');

        end

    end


    %compare to mapping in each condition
    locus_idx=ismember(mapping_input.index,10191);
    
    %match common names
    mapping_common_name=cell(height(mapping_input),1);
    for i=1:height(mapping_input)
        
        temp_str=strsplit(mapping_input.commonName{i},';');
        mapping_common_name{i}=temp_str{1};
        
    end
    
    mapping_beta=nan(length(protein_names),1);
    for i=1:length(protein_names)
    
        protein_idx=ismember(mapping_common_name,protein_names{i});
    
        temp_idx=logical(protein_idx.*locus_idx);
    
        if sum(temp_idx)>0
    
            mapping_beta(i)=mapping_input.beta(temp_idx);
    
        end
    
    end

    
    for i=1:length(condition_names)
    
        v1=strain_volcano_effect{i};
        v2=mapping_beta;
        
        %concordant/non-concordant
        n1=sum((v1<1).*(v2>0));
        n2=sum((v1>1).*(v2>0));
        %RM coherent
        rm_coherent(i)=n2;
        n3=sum((v1<1).*(v2<0));
        %YJM coherent
        yjm_coherent(i)=n3;
        n4=sum((v1>1).*(v2<0));
        
        temp_table=table([n1;n2],[n3;n4],...
            'VariableNames',{'mapping RM up','mapping RM down'},'RowNames',{'ms YJM up','ms YJM down'});
        [h,p,stats]=fishertest(temp_table);
        %text(0.8,-0.4,num2str(p))       
       
        p_coherent(i)=p;
        
    
    
    end

    
    for i=1:length(condition_names)
    
        glc_idx=ismember(v_condition_all,'Glucose');
        condition_idx=ismember(v_condition_all,condition_names{i});
        strain_idx=ismember(v_strain_all,'YJM');
    
        condition_volcano_p{i}=nan(length(protein_names),1);
        condition_volcano_effect{i}=nan(length(protein_names),1);
        
        for j=1:length(protein_names)
    
            temp_idx1=logical(condition_idx.*strain_idx);
            v_temp1=clean_mat(j,temp_idx1);
    
            temp_idx2=logical(glc_idx.*strain_idx);
            v_temp2=clean_mat(j,temp_idx2);
    
            [h,condition_volcano_p{i}(j)]=ttest2(v_temp1,v_temp2);
            condition_volcano_effect{i}(j)=mean(v_temp1,'omitnan')/...
                mean(v_temp2,'omitnan');
    
        end

        v1=log2(condition_volcano_effect{i});
        v2=-log10(mafdr(condition_volcano_p{i},'BHFDR','true'));
    
        n_diff_exp(i)=sum(v2>2);
        
    end
    

    
end


