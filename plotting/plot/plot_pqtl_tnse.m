function [] = plot_pqtl_tnse(gene_name,gene_to_highlight,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names,strain_index]=...
        parse_raw_abundance(dependency_directory,output_directory);

    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);

    %set up F6 haploid matrix
    f6_idx=find(ismember(strain_names,'F6'));
    f6_genotype_idx=strain_index(f6_idx);
    
    [sortedf6_genotype_idx,sort_idx]=sort(f6_genotype_idx,'ascend');
    
    f6_mat=input_mat(:,f6_idx(sort_idx));


    mat_for_tsne=f6_mat;
    
    %zscore
    for i=1:length(orf_names)
        
        v_temp=mat_for_tsne(i,:);
        
        mat_for_tsne(i,:)=(v_temp-mean(v_temp,'omitnan'))./std(v_temp,[],'omitnan');
        
    end

    mat_for_tsne(isnan(mat_for_tsne))=0;
    
    rng(0)
    Y=tsne(mat_for_tsne,'Distance','correlation');

    v_color=zeros(length(orf_names),1);
        
    temp_idx1=logical(ismember(all_pqtl_data.common1,gene_name)+...
        ismember(all_pqtl_data.common2,gene_name));
    
    temp_genes=all_pqtl_data.protein(temp_idx1);
    temp_sign=all_pqtl_data.beta(temp_idx1)>0;
    
    if length(temp_genes)>0
        
        temp_idx2=ismember(orf_names,temp_genes(temp_sign));
        v_color(temp_idx2)=1;
        
        temp_idx2=ismember(orf_names,temp_genes(~temp_sign));
        v_color(temp_idx2)=2;
        
        hold on
        gscatter(Y(:,1),Y(:,2),v_color,[grey;blue;orange],'...',[20 30 30])%'kbr')
        
        plot_idx=ismember(orf_names,gene_to_highlight);
        scatter(Y(plot_idx,1),Y(plot_idx,2),30,'k','filled')
        axis square
        xlim([-20 20])
        ylim([-20 20])
        title(gene_name)
        legend('off')
        axis off
        
    end
        

end


