function [bar_mat] = calculate_biogrid_overlaps(dependency_directory,output_directory)


    %compare trans pQTL connections to various measures [complex and biogrid,
    %etc]

    %trans pQTL pairs
    pqtn_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    trans_idx=pqtn_data.dist>1e3;
    od_idx=pqtn_data.index==0;

    pqtl_to_use=pqtn_data(logical(trans_idx.*~od_idx),:);

    
    
    
    
    
    %biogrid
    biogrid_input=readtable([dependency_directory 'BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.207.tab3.txt']);


    genesA=biogrid_input.SystematicNameInteractorA;
    genesB=biogrid_input.SystematicNameInteractorB;
    v_type=biogrid_input.ExperimentalSystemType;


    query1=pqtl_to_use.protein;

    query2=pqtl_to_use.gene1;
    query3=pqtl_to_use.gene2;

    is_biogrid=zeros(height(pqtl_to_use),1);

    for i=1:length(is_biogrid)

        if mod(i,100)==0
            i
        end

        v1=ismember(genesA,query1{i});
        v2=ismember(genesB,query1{i});

        if sum(v1)>0

            sum1=sum(ismember(genesB(v1),query2{i}));
            sum3=sum(ismember(genesB(v1),query3{i}));

        else

            sum1=0;
            sum3=0;

        end

        if sum(v2)>0

            sum2=sum(ismember(genesB(v2),query2{i}));
            sum4=sum(ismember(genesB(v2),query3{i}));

        else

            sum2=0;
            sum4=0;

        end

        if (sum1+sum2+sum3+sum4)>0

            is_biogrid(i)=1;

        end

    end


    %look up interaction types for hits
    to_lookup=find(is_biogrid);

    for i=1:length(to_lookup)

        query1=pqtl_to_use.protein{to_lookup(i)};

        v1=ismember(genesA,query1);
        v2=ismember(genesB,query1);

        query2=pqtl_to_use.gene1{to_lookup(i)};
        query3=pqtl_to_use.gene2{to_lookup(i)};

        v3=ismember(genesB,query2);
        v4=ismember(genesA,query2);

        v5=ismember(genesB,query3);
        v6=ismember(genesA,query3);

        temp_idx=find((v1.*v3)+(v2.*v4)+(v1.*v5)+(v2.*v6));

        interaction_types{i}=unique(v_type(temp_idx));

    end


    for i=1:length(interaction_types)

        if length(interaction_types{i})==1

            is_genetic(i)=strcmp(interaction_types{i},'genetic');
            is_physical(i)=strcmp(interaction_types{i},'physical');

        elseif length(interaction_types{i})==2

            is_both(i)=1;

        end

    end



    clear bar_mat
    %biogrid; then split by biogrid type
    %bar_mat(1,1)=sum(is_complex)/length(is_complex);
    bar_mat(1,1)=sum(is_biogrid)/length(is_biogrid);
    bar_mat(1,2)=1-bar_mat(1,1);
    bar_mat(1,3)=0;

    bar_mat(3,1)=sum(is_genetic)/length(is_genetic);
    bar_mat(3,2)=sum(is_physical)/length(is_genetic);
    bar_mat(3,3)=sum(is_both)/length(is_genetic);
    bar_mat(3,4)=1-bar_mat(3,3)-bar_mat(3,2)-bar_mat(3,1);


    



end