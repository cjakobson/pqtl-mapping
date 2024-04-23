%plot effect of pQTL
function []=plot_locus_effect(gene_name,common_name,locus_number,y_lim,...
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

    to_plot{1}=v_temp(rm_idx);
    to_plot{2}=v_temp(yjm_idx);
    %to_plot{3}=v_temp(sgrp_idx);
    
    easy_box([to_plot f6_for_box_plot]);
    ylim([0 y_lim])
    title(common_name)
    
    
end


