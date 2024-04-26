function [rm_od,yjm_od,strain_od,f6_od,...
    rm_idx,yjm_idx,strain_idx]=parse_od(dependency_directory,output_directory)

    %data for heritability calcs
    [num, txt]=xlsread([dependency_directory '50-0040_segregants_meta-table_v6.xls']);

    for i=1:length(txt)

        name_to_parse=txt{i,1};

        temp_str=strsplit(name_to_parse,{'_',','});

        if strcmp(temp_str{1},'F6')

            strain_idx(i)=str2num(temp_str{3});
            strain_od(i)=num(i-1);

            %output parental data for heritability analysis
        elseif strcmp(temp_str{1},'RM11')

            rm_idx{i}=[temp_str{3:5}];
            rm_od(i)=num(i-1);

        elseif strcmp(temp_str{1},'YJM975')

            yjm_idx{i}=[temp_str{3:5}];
            yjm_od(i)=num(i-1);

        else

            strain_idx(i)=nan;
            strain_od(i)=nan;

        end

    end

    strain_od(isnan(strain_idx))=[];
    strain_idx(isnan(strain_idx))=[];

    strain_od(strain_idx==0)=[];
    strain_idx(strain_idx==0)=[];

    f6_od(strain_idx)=strain_od;

    rm_od(cellfun(@isempty,rm_idx))=[];
    rm_idx(cellfun(@isempty,rm_idx))=[];

    yjm_od(cellfun(@isempty,yjm_idx))=[];
    yjm_idx(cellfun(@isempty,yjm_idx))=[];


end

