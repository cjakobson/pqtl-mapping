function [h_squared_rm,h_squared_yjm,h_squared_mean,var_exp]=...
        calculate_heritability(dependency_directory,output_directory)
    
    
    mapping_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv'],...
        'FileType','text');
    
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);
    
    
    [rm_od,yjm_od,~,~,...
        rm_idx,yjm_idx,~]=...
        parse_od(dependency_directory,output_directory);
    
    
    for i=1:length(orf_names)

        clear v_od
        clear v_protein
        for j=1:length(rm_idx)

            temp_idx_2=find(ismember(strain_merge_idx,rm_idx{j}));
            v_od(j)=rm_od(j);
            v_protein(j)=input_mat(i,temp_idx_2);

        end

        if sum(~isnan(v_protein)>0)
            temp_tbl=table(v_protein',v_od');
            temp_lme=fitlme(temp_tbl,'Var1~Var2');
            h_squared_rm(i)=1-sqrt(temp_lme.MSE)/mean(v_protein,'omitnan');
        end

        clear v_od
        clear v_protein
        for j=1:length(yjm_idx)

            temp_idx_2=find(ismember(strain_merge_idx,yjm_idx{j}));

            if ~isempty(temp_idx_2)
                v_od(j)=yjm_od(j);
                v_protein(j)=input_mat(i,temp_idx_2);
            end

        end

        if sum(~isnan(v_protein))>2
            temp_tbl=table(v_protein',v_od');
            temp_lme=fitlme(temp_tbl,'Var1~Var2');
            h_squared_yjm(i)=1-sqrt(temp_lme.MSE)/mean(v_protein,'omitnan');
        end


        orfIdx=ismember(mapping_data.protein,orf_names{i});
        if sum(orfIdx)>0
            var_exp(i)=sum(mapping_data.varExp(orfIdx));
        end
        
    end


    h_squared_mean=mean([h_squared_rm; h_squared_yjm]);
    
    
    
end


