function [prot_corr, prot_corr_1K, prot_corr_5K] =...
    calculate_correlations(dependency_directory,output_directory)


    %protein-protein correlations

    %1K proteomes data
    input_data1K=readtable([dependency_directory '230502_proteomes.csv'],'FileType','text');

    orf_names_1K=input_data1K.Protein_Group;
    abundance_mat_1K=table2array(input_data1K(:,2:(end-1)));

    %5K proteomes data
    input_data5K=readtable([dependency_directory 'yeast5k_impute_wide.csv'],...
        'FileType','text');

    orf_names_5K_pdb=input_data5K.Protein_Group;
    abundance_mat_5K=table2array(input_data5K(:,2:end));

    %convert 5K ORFs to systematic names

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


    for i=1:length(orf_names_5K_pdb)

        temp_idx=ismember(pdb_id,orf_names_5K_pdb{i});

        if sum(temp_idx)>0

            orf_names_5K{i}=systematic_name_pdb{temp_idx};

        else

            orf_names_5K{i}='NA';

        end

    end

    orf_names_5K=orf_names_5K';



    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);

    f6_mat=input_mat(:,f6_idx);


    %assess protein-protein correlations
    prot_corr=nan(length(orf_names));
    prot_corr_1K=nan(length(orf_names));
    prot_corr_5K=nan(length(orf_names));
    m=1;
    n=1;
    o=1;
    for i=1:length(orf_names)


        v1=f6_mat(i,:)';

        for j=(i+1):length(orf_names)

            v2=f6_mat(j,:)';
            [r p]=corr(v1,v2,'rows','complete');

            prot_corr(i,j)=r;

        end


        %do 1K
        idx_1K=find(ismember(orf_names_1K,orf_names{i}));

        if ~isempty(idx_1K)

            v3=abundance_mat_1K(idx_1K,:)';

            for j=(i+1):length(orf_names)

                idx_1K=find(ismember(orf_names_1K,orf_names{j}));

                if ~isempty(idx_1K)

                    v4=abundance_mat_1K(idx_1K,:)';

                    [r p]=corr(v3,v4,'rows','complete');

                    prot_corr_1K(i,j)=r;

                end

            end

        end

        %do 5K
        idx_5K=find(ismember(orf_names_5K,orf_names{i}));

        if ~isempty(idx_5K)

            v5=abundance_mat_5K(idx_5K,:)';

            for j=(i+1):length(orf_names)

                idx_5K=find(ismember(orf_names_5K,orf_names{j}));

                if ~isempty(idx_5K)

                    v6=abundance_mat_5K(idx_5K,:)';

                    [r p]=corr(v5,v6,'rows','complete');

                    prot_corr_5K(i,j)=r;

                end

            end

        end

    end



end