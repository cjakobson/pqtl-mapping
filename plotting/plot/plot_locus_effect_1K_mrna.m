function []=plot_locus_effect_1K_mrna(gene_name,common_name,locus_number,y_lim,...
    dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    
    %genotype data
    load([dependency_directory 'phasedGenotype.mat'])
    
    %proteome data
    protein_data=readtable([dependency_directory '211031_SegregantProteomicsData_DetectionThreshold80_genes_ORF.tsv'],...
        'FileType','text');

    protein_mat=table2array(protein_data(:,2:end));

    orf_names=table2array(protein_data(:,1));

    temp_strain=protein_data.Properties.VariableNames(2:end);
    m=1;
    n=1;
    o=1;
    for i=1:length(temp_strain)

        temp_str=strsplit(temp_strain{i},'_');

        if strcmp(temp_str{1},'F6')

            f6_idx(i)=str2num(temp_str{3});

        else

            f6_idx(i)=nan;

        end

        if strcmp(temp_str{1},'RM11')

            rm_idx(m)=i;
            m=m+1;

        elseif strcmp(temp_str{1},'YJM975')

            yjm_idx(n)=i;
            n=n+1;

        elseif ~strcmp(temp_str{1},'F6')

            sgrp_idx(o)=i;
            o=o+1;

        end

    end

    
    qtl_to_use=locus_number;
    
    clear qtl_genotype
    qtl_genotype{1,1}=find(phasedGenotype(:,qtl_to_use)==1);
    qtl_genotype{2,1}=find(phasedGenotype(:,qtl_to_use)==-1);


    clear genotypes_for_box_plot
    m=1;
    for j=1:2
        v1=qtl_genotype{j,1};
        genotypes_for_box_plot{m}=v1;
        m=m+1;
    end


    v_temp=protein_mat(ismember(orf_names,gene_name),:);

    clear f6_for_box_plot
    for j=1:2
        f6_for_box_plot{j}=v_temp(ismember(f6_idx,genotypes_for_box_plot{j}));
    end

%     to_plot{1}=v_temp(rm_idx);
%     to_plot{2}=v_temp(yjm_idx);
    
    to_plot{1}=f6_for_box_plot{1};
    to_plot{2}=f6_for_box_plot{2};

    %normalize to RM median
    temp_median=median(to_plot{1},'omitnan');
    
    to_plot{1}=to_plot{1}/temp_median;
    to_plot{2}=to_plot{2}/temp_median;
    
    
    cis_pqtn_data=calculate_1K_replication_mrna(dependency_directory,output_directory);
    
    temp_idx=find(ismember(cis_pqtn_data.protein,gene_name));
    to_plot{3}=cis_pqtn_data.rm_array{temp_idx};
    to_plot{4}=cis_pqtn_data.yjm_array{temp_idx};
    
    temp_median=median(to_plot{3},'omitnan');
    
    to_plot{3}=to_plot{3}/temp_median;
    to_plot{4}=to_plot{4}/temp_median;
    
    
    hold on
    easy_box(to_plot)
    title([cis_pqtn_data.commonName{temp_idx} ' mRNA'])
    xticks(1:4)
    xtickangle(45)
    xticklabels({'RM allele','YJM allele','RM allele','YJM allele'})
    ylim([0 y_lim])
    v_temp=ylim;
    ylim([0 v_temp(2)])
    text(1,0.8*v_temp(2),['\beta = ' num2str(cis_pqtn_data.beta(temp_idx))])
    for j=3:4
        text(j,0.1*v_temp(2),num2str(length(to_plot{j})))
    end
    
    
end


