function [pqtn_sorted, qtn_sorted] = ...
    calculate_pqtn_qtn_scores(gene_name,locus_to_use,dependency_directory,output_directory)

    %plot pQTN scores for all its targets
    pqtl_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    pqtl_idx_to_use=pqtl_input.protein(logical(ismember(pqtl_input.gene1,gene_name)+...
        ismember(pqtl_input.gene2,gene_name)));

    pqtl_loci_to_use=(locus_to_use-10):(locus_to_use+10);
    pqtl_loci_to_use=pqtl_loci_to_use+1;  %OD600 variable offset


    m=1;
    clear v_to_plot v_max
    for i=1:length(pqtl_idx_to_use)

        load([dependency_directory 'linearPqtl/' pqtl_idx_to_use{i} '.mat'])

        for j=1:length(pqtl_loci_to_use)

            if sum(ismember(posToMap,pqtl_loci_to_use(j)))>0

                idx_to_use=find(ismember(posToMap,pqtl_loci_to_use(j)));

                temp_mat=ph2{idx_to_use};
                %subplot(3,5,j)
                %hold on
                temp_mat(temp_mat==-1)=NaN;
                v_to_plot{m}=min(real(-log10(temp_mat)),[],2);
                v_max(m)=max(v_to_plot{m});
                m=m+1;

            end

        end

    end

    %sort from back to front
    [~,sort_idx]=sort(v_max,'descend');
    
    pqtn_sorted=v_to_plot(sort_idx);
    
    
    
    
    growth_input=readtable([dependency_directory 'linearNoRad.csv']);

    growth_idx_to_use=logical(ismember(growth_input.gene1,gene_name)+...
        ismember(growth_input.gene2,gene_name));

    temp_table=growth_input(growth_idx_to_use,:);

    %remove hets
    temp_table(temp_table.bPos>(12054+41),:)=[];

    clear v_to_plot v_max
    for n=1:height(temp_table)

        load([dependency_directory 'linearNoRad/' temp_table.time{n} ' '...
            temp_table.condition{n} '-rad.mat'])

        temp_locus=temp_table.bPos(n);

        temp_idx=find(ismember(posToMap,temp_locus));

        temp_mat=ph2{temp_idx};

        temp_mat(temp_mat==-1)=NaN;
        v_to_plot{n}=min(real(-log10(temp_mat)),[],2);
        v_max(n)=max(v_to_plot{n});

    end


    %sort from back to front
    [~,sort_idx]=sort(v_max,'descend');
    
    qtn_sorted=v_to_plot(sort_idx);

    

end