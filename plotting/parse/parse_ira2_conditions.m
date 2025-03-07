%parse strain names and raw abundance data
function [mean_mat,v_label_mean,v_strain_mean,v_condition_mean,...
    clean_mat,v_label_all,v_strain_all,v_condition_all,...
    z_mat,protein_names]=parse_ira2_conditions(dependency_directory,output_directory)

    input_data=readtable([dependency_directory '50-0040-diff-condition-proteomes-pass2-all-samples.xlsx']);
    
    input_data(ismember(input_data.strain,''),:)=[];


    %remove extraneous conditions
    input_data(ismember(input_data.condition,'EtOH-harvest-high'),:)=[];
    input_data(ismember(input_data.condition,'Rapamycin-500nM'),:)=[];
    
    
    
    strain_names=unique(input_data.strain);
    
    replicate_names=unique(input_data.replicate);
    
    condition_names=unique(input_data.condition);
    
    
    protein_names=unique(input_data.COMMON);
    protein_names(cellfun(@isempty,protein_names))=[];
    protein_names(ismember(protein_names,'NA'))=[];
    
    
    %include replicates for PCA
    all_mat=nan(length(protein_names),length(strain_names)*length(condition_names)...
        *length(replicate_names));
    
    for i=1:length(protein_names)
    
        temp_table=input_data(ismember(input_data.COMMON,protein_names{i}),:);
    
        m=1;
        for j=1:length(strain_names)
            
            temp_strain_idx=ismember(temp_table.strain,strain_names{j});
    
            for k=1:length(condition_names)
    
                temp_condition_idx=ismember(temp_table.condition,condition_names{k});
    
                for l=1:length(replicate_names)
    
                    temp_replicate_idx=ismember(temp_table.replicate,replicate_names{l});
    
                    temp_idx=logical(temp_strain_idx.*temp_condition_idx.*temp_replicate_idx);
    
                    if sum(temp_idx)>0
        
                        all_mat(i,m)=temp_table.abundance(temp_idx);
        
                    end
        
                    v_label_all{m}=[strain_names{j} '_' condition_names{k}...
                        '_' replicate_names{l}];
    
                    v_strain_all{m}=strain_names{j};
                    v_condition_all{m}=condition_names{k};
        
                    m=m+1;
    
                end
    
            end
    
        end
    
    end
    
    
    outlier_idx=[4,44,46,120];

    clean_mat=all_mat;
    clean_mat(:,outlier_idx)=nan;
    
    %repeat PCA
    
    
    %pca
    z_mat=nan(size(clean_mat));
    for i=1:length(protein_names)
    
        v_temp=clean_mat(i,:);
    
        z_mat(i,:)=(v_temp-mean(v_temp,'omitnan'))./std(v_temp,[],'omitnan');
    
    end
    
    z_mat(isnan(z_mat))=0;



    mean_mat=nan(length(protein_names),length(strain_names)*length(condition_names));
    %mean for downstream analyses
    for i=1:length(protein_names)
    
        m=1;
        for j=1:length(strain_names)
            
            temp_strain_idx=ismember(v_strain_all,strain_names{j});
    
            for k=1:length(condition_names)
    
                temp_condition_idx=ismember(v_condition_all,condition_names{k});
    
                temp_idx=logical(temp_strain_idx.*temp_condition_idx);

                mean_mat(i,m)=mean(clean_mat(i,temp_idx),'omitnan');
                
    
                v_label_mean{m}=[strain_names{j} '_' condition_names{k}];

                v_strain_mean{m}=strain_names{j};
                v_condition_mean{m}=condition_names{k};
    
                m=m+1;

                
    
            end
    
        end
    
    end
    

end

