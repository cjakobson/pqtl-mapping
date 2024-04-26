function [physical_gene_order,breaks_to_plot,breaks_to_plot2]=...
    calculate_chr_breaks(dependency_directory,output_directory)

    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);


    all_pqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    all_pqtl(all_pqtl.bPos==1,:)=[];

    orf_names=unique(all_pqtl.protein);

    [sorted_orf_names,alpha_idx]=sort(orf_names);


    %find chromosome breaks
    m=0;
    for i=2:length(sorted_orf_names)

        break_idx(i)=sorted_orf_names{i}(2)~=sorted_orf_names{i-1}(2);
        if break_idx(i)
            m=m+1;
        end
        chr_idx(i)=m;

    end

    chr_breaks=find(break_idx);

    chr_names={'M','I','II','III','IV','V','VI','VII','VIII','IX','X',...
        'XI','XII','XIII','XIV','XV','XVI'};

    physical_gene_order=[];

    %separate L and R arm and reorder
    for i=0:17

        temp_subset=sorted_orf_names(chr_idx==i);

        %reverse order of left arm
        is_left=zeros(length(temp_subset),1);
        for j=1:length(temp_subset)

            if length(temp_subset{j})>=3
                is_left(j)=temp_subset{j}(3)=='L';
            end

        end

        if sum(is_left)>0

            temp_subset(1:sum(is_left))=flipud(temp_subset(1:sum(is_left)));

        end

        physical_gene_order=[physical_gene_order; temp_subset];

    end

    %to annotate chrs on plots
    %variant_info=readtable('variant_info.csv');
    v_diff=variant_info.chr(2:end)~=variant_info.chr(1:(end-1));
    breaks_to_plot=find(v_diff);

    %in other dimension
    breaks_to_plot2=find(break_idx);
    breaks_to_plot2=breaks_to_plot2(2:end);

end
