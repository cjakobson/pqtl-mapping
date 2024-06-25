function []=plot_correlation_heatmap(gene_list,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,orf_names,strain_index]=...
            parse_raw_abundance(dependency_directory,output_directory);

    mat_to_use1=input_mat(:,rm_idx);
    mat_to_use2=input_mat(:,yjm_idx);
        

    %filter missing proteins
    genes_to_use=gene_list(ismember(gene_list,orf_names));

    corr_mat=zeros(length(gene_list));
    for i=1:length(gene_list)

        temp_idx1=ismember(orf_names,gene_list{i});

        for j=1:length(gene_list)

            temp_idx2=ismember(orf_names,gene_list{j});

            v1=mat_to_use1(temp_idx1,:)';
            v2=mat_to_use1(temp_idx2,:)';
            [r p]=corr(v1,v2,'rows','complete');

            corr_mat1(i,j)=r;
            
            v1=mat_to_use2(temp_idx1,:)';
            v2=mat_to_use2(temp_idx2,:)';
            [r p]=corr(v1,v2,'rows','complete');

            corr_mat2(i,j)=r;

        end

    end
    
    
    imagesc((corr_mat1+corr_mat2)/2,[-0.75 0.75])
    axis square
    colormap gray
    colorbar
    yticklabels(gene_list)
    xticks(0.5:1:length(gene_list))
    yticks(xticks)
    grid on

    m = size(get(gcf,'colormap'),1);
    %red to blue colormap
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    c = [r g b]; 
    %colormap(flipud(c))
    colormap(c)



end


