function [gene_name,gene_chr,gene_start,gene_end]=get_gene_coords(dependency_directory)

fasta_input=fastaread([dependency_directory 'orf_genomic_all.fasta']);

chr_array={'I','II','III','IV','V','VI','VII','VIII','IX','X',...
    'XI','XII','XIII','XIV','XV','XVI'};
for i=1:length(fasta_input)

    temp_str=strsplit(fasta_input(i).Header,',');

    temp_str1=strsplit(temp_str{1},' ');
    gene_name{i}=temp_str1{1};

    temp_str2=strsplit(temp_str{2},' ');
    temp_chr=temp_str2{3};

    temp_chr_num=find(ismember(chr_array,temp_chr));

    if ~isempty(temp_chr_num)
        gene_chr(i)=temp_chr_num;
    else    %mito
        gene_chr(i)=17;
    end

    temp_pos_str=temp_str2{5};
    temp_str3=strsplit(temp_pos_str,'-');

    gene_start(i)=str2num(temp_str3{1});
    gene_end(i)=str2num(temp_str3{2});

end


end


