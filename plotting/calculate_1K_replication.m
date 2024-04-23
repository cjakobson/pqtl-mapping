function cis_pqtn_data=calculate_1K_replication(dependency_directory,output_directory)

    %all segregating
    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

    
    %1002 genome and proteome data
    proteomics_1002_data=readtable([dependency_directory '211025_ProteomicsData_SMmedia_DIA-NN-1.8_genes_ORF.tsv'],...
        'FileType','text');

    load([dependency_directory '1002data.mat']);

    
    %parse strain names--some are missing from proteomic data
    v_proteomic_strains=proteomics_1002_data.Properties.VariableNames(2:end);
    for i=1:length(v_proteomic_strains)
        if strcmp(v_proteomic_strains{i},'Reference')
            proteomic_genotype1(:,i)=zeros(height(variant_info),1);
            proteomic_genotype2(:,i)=zeros(height(variant_info),1);
        elseif sum(ismember(strainString,v_proteomic_strains{i}))>0
            proteomic_genotype1(:,i)=minGenotype(:,ismember(strainString,v_proteomic_strains{i}));
            proteomic_genotype2(:,i)=minGenotype2(:,ismember(strainString,v_proteomic_strains{i}));
        end
    end

    %RM11 and YJM975 codes
    rm_name='AAA';
    rm_idx=find(ismember(v_proteomic_strains,rm_name));

    yjm_name='SACE_YBA';
    yjm_idx=find(ismember(v_proteomic_strains,yjm_name));

    
    all_pqtn_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    dist_thresh=1e3;
    cis_pqtn_data=all_pqtn_data(all_pqtn_data.dist<dist_thresh,:);
    
    %plot effects of cis acting variants
    for i=1:height(cis_pqtn_data)

        protein_query=cis_pqtn_data.protein{i};

        protein_idx=ismember(proteomics_1002_data.Protein_Group,protein_query);

        if sum(protein_idx)>0

            qtl_idx=cis_pqtn_data.index(i);

            temp_genotype1=proteomic_genotype1(qtl_idx,:);
            temp_genotype2=proteomic_genotype2(qtl_idx,:);

            %code as homozygous RM/het/homozygous YJM only to start
            v_genotype=nan(1,length(v_proteomic_strains));

            v_rm=find(logical((temp_genotype1==temp_genotype1(rm_idx)).*...
                (temp_genotype2==temp_genotype2(rm_idx))));
            %v_rm(v_rm==rm_idx)=[]; %remove RM itself
            v_genotype(v_rm)=1;

            v_yjm=find(logical((temp_genotype1==temp_genotype1(yjm_idx)).*...
                (temp_genotype2==temp_genotype2(yjm_idx))));
            %v_yjm(v_yjm==yjm_idx)=[];
            v_genotype(v_yjm)=-1;

            v_het=find(logical(((temp_genotype1==temp_genotype1(rm_idx)).*...
                (temp_genotype2==temp_genotype2(yjm_idx)))+...
                ((temp_genotype1==temp_genotype1(yjm_idx)).*...
                (temp_genotype2==temp_genotype2(rm_idx)))));
            v_genotype(v_het)=0;

            cis_pqtn_data.rm_array{i}=table2array(proteomics_1002_data(protein_idx,v_rm+1));
            cis_pqtn_data.het_array{i}=table2array(proteomics_1002_data(protein_idx,v_het+1));
            cis_pqtn_data.yjm_array{i}=table2array(proteomics_1002_data(protein_idx,v_yjm+1));

            cis_pqtn_data.rm_mean(i)=mean(cis_pqtn_data.rm_array{i});
            cis_pqtn_data.yjm_mean(i)=mean(cis_pqtn_data.yjm_array{i});
            cis_pqtn_data.het_mean(i)=mean(cis_pqtn_data.het_array{i});

            cis_pqtn_data.rm_freq(i)=length(cis_pqtn_data.rm_array{i});
            cis_pqtn_data.yjm_freq(i)=length(cis_pqtn_data.yjm_array{i});
            cis_pqtn_data.het_freq(i)=length(cis_pqtn_data.het_array{i});

            cis_pqtn_data.maf(i)=min([cis_pqtn_data.rm_freq(i) cis_pqtn_data.yjm_freq(i) cis_pqtn_data.het_freq(i)]);

            cis_pqtn_data.rm_rarer(i)=cis_pqtn_data.rm_freq(i)<cis_pqtn_data.yjm_freq(i);
            cis_pqtn_data.yjm_rarer(i)=cis_pqtn_data.yjm_freq(i)<cis_pqtn_data.rm_freq(i);

        end

    end


    cis_pqtn_data(cis_pqtn_data.maf<0,:)=[];
    
end
