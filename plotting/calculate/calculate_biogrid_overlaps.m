function [bar_mat] = calculate_biogrid_overlaps(dependency_directory,output_directory)


    %compare trans pQTL connections to various measures [complex and biogrid,
    %etc]

    %trans pQTL pairs
    pqtn_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    trans_idx=pqtn_data.dist>1e3;
    od_idx=pqtn_data.index==0;

    pqtl_to_use=pqtn_data(logical(trans_idx.*~od_idx),:);

    
    %use pre-calculated biogrid matrices
    load([dependency_directory 'biogrid_data.mat'])  
    
    is_biogrid=zeros(height(pqtl_to_use),1);

    for i=1:length(is_biogrid)

        target_idx=find(ismember(all_genes,pqtl_to_use.protein{i}));

        gene_idx=find(ismember(all_genes,{pqtl_to_use.gene1{i},pqtl_to_use.gene2{i}}));

        if ~isempty(target_idx)&&~isempty(gene_idx)

            is_biogrid(i)=logical(sum(interaction_mat(target_idx,gene_idx)));

        end

    end


    load([dependency_directory 'biogrid_data_genetic.mat'])  
    
    is_genetic=zeros(height(pqtl_to_use),1);

    for i=1:length(is_genetic)

        target_idx=find(ismember(all_genes,pqtl_to_use.protein{i}));

        gene_idx=find(ismember(all_genes,{pqtl_to_use.gene1{i},pqtl_to_use.gene2{i}}));

        if ~isempty(target_idx)&&~isempty(gene_idx)

            is_genetic(i)=logical(sum(interaction_mat_genetic(target_idx,gene_idx)));

        end

    end

    load([dependency_directory 'biogrid_data_physical.mat'])  
    
    is_physical=zeros(height(pqtl_to_use),1);

    for i=1:length(is_genetic)

        target_idx=find(ismember(all_genes,pqtl_to_use.protein{i}));

        gene_idx=find(ismember(all_genes,{pqtl_to_use.gene1{i},pqtl_to_use.gene2{i}}));

        if ~isempty(target_idx)&&~isempty(gene_idx)

            is_physical(i)=logical(sum(interaction_mat_physical(target_idx,gene_idx)));

        end

    end

    is_both=logical(is_physical.*is_genetic);

    %don't double count
    is_physical(is_both)=0;
    is_genetic(is_both)=0;


    sum(is_biogrid)
    sum(is_genetic)/sum(is_biogrid)
    sum(is_physical)/sum(is_biogrid)
    sum(is_both)/sum(is_biogrid)


    clear bar_mat
    %biogrid; then split by biogrid type
    %bar_mat(1,1)=sum(is_complex)/length(is_complex);
    bar_mat(1,1)=sum(is_biogrid)/length(is_biogrid)
    bar_mat(1,2)=1-bar_mat(1,1);
    bar_mat(1,3)=0;

    bar_mat(3,1)=sum(is_genetic)/sum(is_biogrid);
    bar_mat(3,2)=sum(is_physical)/sum(is_biogrid);
    bar_mat(3,3)=sum(is_both)/sum(is_biogrid);
    bar_mat(3,4)=1-bar_mat(3,3)-bar_mat(3,2)-bar_mat(3,1);


    sum(sum(triu(interaction_mat,1)))
    sum(sum(triu(interaction_mat_physical,1)))/sum(sum(triu(interaction_mat,1)))
    sum(sum(triu(interaction_mat_genetic,1)))/sum(sum(triu(interaction_mat,1)))

    %binomial test vs all of biogrid
    p=binocdf(sum(is_genetic)+sum(is_both),sum(is_biogrid),...
        sum(sum(triu(interaction_mat_genetic,1)))/sum(sum(triu(interaction_mat,1))))

    p=binocdf(sum(is_physical)+sum(is_both),sum(is_biogrid),...
        sum(sum(triu(interaction_mat_physical,1)))/sum(sum(triu(interaction_mat,1))))


end