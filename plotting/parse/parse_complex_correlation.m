function [corr_string_input_data,complex_corr]=...
    parse_complex_correlation(dependency_directory,output_directory)

    corr_string_input_data=readtable([dependency_directory 'corrData_withString_withCellmap.csv']);

    %plot correlations within/outside complexes
    complex_data=readtable([dependency_directory '20230508_ComplexEBI_saccharomyces_cerevisiae_slim_JH.csv'],...
        'FileType','text');

    pdb_dictionary=tdfread([dependency_directory 'uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.07.07-20.59.31.22.tsv']);

    pdb_id=cell(length(pdb_dictionary.Entry),1);
    systematic_name_pdb=pdb_id;
    m=1;
    for i=1:length(pdb_dictionary.Entry)

        temp_str=pdb_dictionary.Entry(i,:);
        temp_str(temp_str==' ')=[];
        pdb_id{m}=temp_str;

        temp_str=pdb_dictionary.Gene_Names_0x28ordered_locus0x29(i,:);
        temp_str(temp_str==' ')=[];
        systematic_name_pdb{m}=temp_str;
        m=m+1;

    end


    clear complex_corr
    clear systematic_query
    n=1;
    for i=1:height(complex_data)

        temp_str=strsplit(complex_data.proteins{i},'|');

        %look up uniprot IDs
        m=1;
        for j=1:length(temp_str)

            temp_idx=find(ismember(pdb_id,temp_str{j}));

            if ~isempty(temp_idx)

                systematic_query{i}{m}=systematic_name_pdb{temp_idx};
                m=m+1;

            end

        end

        if m>1

            for j=1:length(systematic_query{i})

                temp_idx1=ismember(corr_string_input_data.v1,systematic_query{i}{j});

                if sum(temp_idx1)>0

                    for k=1:length(systematic_query{i})

                        temp_idx2=ismember(corr_string_input_data.v2,systematic_query{i}{k});

                        temp_idx=logical(temp_idx1.*temp_idx2);

                        if sum(temp_idx)>0

                            complex_corr(n)=corr_string_input_data.v3(temp_idx);
                            n=n+1;

                        end

                    end

                end

            end

        end

    end
    
    
    

end



