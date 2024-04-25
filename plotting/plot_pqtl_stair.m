function []= plot_pqtl_stair(gene_name,pqtls,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];

    load([dependency_directory 'phasedGenotype.mat'])

    for i=1:length(pqtls)
    
        qtl_genotype{1,i}=find(phasedGenotype(:,pqtls(i))==-1);
        qtl_genotype{2,i}=find(phasedGenotype(:,pqtls(i))==1);
    
    end
    
    m=1;
    for i=1:2
        v1=qtl_genotype{i,1};
        for j=1:2
            v2=qtl_genotype{j,2};
            for k=1:2
                v3=qtl_genotype{k,3};
                for l=1:2
                    v4=qtl_genotype{l,4};
                    v_temp=intersect(v1,v2);
                    v_temp=intersect(v_temp,v3);
                    v_temp=intersect(v_temp,v4);
                    genotypes_for_boxplot{m}=v_temp;
                    m=m+1;
                end
            end
        end
    end
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);

    v_temp=input_mat(ismember(orf_names,gene_name),:);

    for i=1:length(genotypes_for_boxplot)
        f6_for_boxplot{i}=v_temp(ismember(strain_index,genotypes_for_boxplot{i}));
    end
    
    to_plot{1}=v_temp(rm_idx);
    to_plot{2}=v_temp(yjm_idx);
    %toPlot{3}=vTemp(sgrpIdx);
    
    %norm to YJM median
    for i=1:length(f6_for_boxplot)
        f6_for_boxplot{i}=f6_for_boxplot{i}./median(to_plot{2},'omitnan');
    end
    for i=1:length(to_plot)
        to_plot{i}=to_plot{i}./median(to_plot{2},'omitnan');
    end
    
    hold on
    easy_box([f6_for_boxplot to_plot]);
    ylim([0 4.5])
    title('Mcr1')



end


