function [] = plot_niche_example(chr_to_plot,pos_to_plot,gene_label,...
    dependency_directory,output_directory)

    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    niche_input=readtable([dependency_directory '1K_common_annotated_niche_all.csv']);
    
    niche_names={'Baking','Dairy','Ferm','Human','Other','Plant','Soil','Waste'};
    
    idx_to_plot=find(ismember(niche_input.chr,chr_to_plot).*niche_input.pos==pos_to_plot);
    
    v1=table2array(niche_input(idx_to_plot,9:16)); %ref
    v2=table2array(niche_input(idx_to_plot,17:24)); %   alt
    
    v3=v2./(v1+v2);
    
    hold on
    scatter(1:length(v3),v3,50,'k','filled')
    ylim([0 1])
    xlim([0.5 length(niche_names)+0.5])
    xticks(1:length(niche_names))
    xtickangle(45)
    xticklabels(niche_names)
    title(gene_label)
    ylabel('alt allele freq.')
    text(niche_input.niche_id(idx_to_plot),0.5,num2str(niche_input.niche_q_value(idx_to_plot)))
    
    %change this to all other, not overall
    to_use=1:length(v1);
    to_use(niche_input.niche_id(idx_to_plot))=[];
    overall_alt=sum(v2(to_use))/(sum(v1(to_use))+sum(v2(to_use)));
    plot(xlim,[overall_alt overall_alt],':r')

end


