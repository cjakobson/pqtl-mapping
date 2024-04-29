function [] = plot_pqtl_5K(gene_name,dependency_directory,output_directory)

    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    idx_to_use=find(ismember(all_pqtl_data.gene1,gene_name)+...
        ismember(all_pqtl_data.gene2,gene_name));
    
    
    %5K data
    deletion_data=readtable([dependency_directory 'yeast5k_impute_wide.csv']);

    %extract gene names
    deletion_mat=table2array(deletion_data(:,2:end));


    deletion_uniprot=deletion_data.Protein_Group;


    %zscore proteins across all deletions
    for i=1:length(deletion_uniprot)

        v_temp=deletion_mat(i,:);
        deletion_z_mat(i,:)=(v_temp-mean(v_temp,'omitnan'))/std(v_temp,[],'omitnan');

    end




    pdb_dict=tdfread([dependency_directory 'uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.07.07-20.59.31.22.tsv']);

    pdb_id=cell(length(pdb_dict.Entry),1);
    systematic_name_pdb=pdb_id;
    m=1;
    for i=1:length(pdb_dict.Entry)

        temp_str=pdb_dict.Entry(i,:);
        temp_str(temp_str==' ')=[];
        pdb_id{m}=temp_str;

        temp_str=pdb_dict.Gene_Names_0x28ordered_locus0x29(i,:);
        temp_str(temp_str==' ')=[];
        systematic_name_pdb{m}=temp_str;
        m=m+1;

    end


    for i=1:length(deletion_uniprot)

        temp_idx=ismember(pdb_id,deletion_uniprot{i});

        if sum(temp_idx)>0

            deletion_protein_systematic{i}=systematic_name_pdb{temp_idx};

        else

            deletion_protein_systematic{i}='NA';

        end

    end

    temp_names=deletion_data.Properties.VariableNames(2:end);
    for i=1:length(temp_names)

        temp_str=strsplit(temp_names{i},'_');
        deletion_gene_systematic{i}=temp_str{5};

    end

    
    clear mapping_beta deletion_z_score
    for i=1:length(idx_to_use)
    
        temp_locus=gene_name;
        temp_target=all_pqtl_data.protein{idx_to_use(i)};
    
        mapping_beta(i)=all_pqtl_data.beta(idx_to_use(i));
        
        temp_idx1=ismember(deletion_gene_systematic,temp_locus);
        temp_idx2=ismember(deletion_protein_systematic,temp_target);
        
        if (sum(temp_idx1)>0)&&(sum(temp_idx2)>0)
    
            v_temp=deletion_z_mat(temp_idx2,temp_idx1);
            deletion_z_score(i)=v_temp(end);
            
        else
            
            deletion_z_score(i)=nan;
    
        end
    
    end
    
    
    
    hold on
    scatter(deletion_z_score,mapping_beta,10,'k','filled')
    axis square
    xlim([-8 8])
    ylim([-0.8 0.8])
    xlabel('proteomics Z score')
    ylabel('mapping \beta')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')
    title(gene_name)
    
    [r p]=corr(deletion_z_score',mapping_beta','rows','complete');
    text(3,0.4,num2str(r))
    text(3,0.3,num2str(p))
    
    
    
end