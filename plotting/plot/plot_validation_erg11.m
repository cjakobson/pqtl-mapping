function [] = plot_validation_erg11(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;

    input_data=readtable([dependency_directory 'validation1_remeasured_metadata-joined_to-chris_v1.csv']);
    
    protein_names=unique(input_data.GENENAME);
    protein_names(cellfun(@isempty,protein_names))=[];

    input_data(ismember(input_data.condition,'flc'),:)=[];


    genes_to_plot={'ERG11'};
    strains_to_plot={'YJM975',...
        'YJM975_ERG11_122014T-C',...
        'YJM975_Erg11_433Asn-Lys',...
        'YJM975_ERG11_122014T-C-Erg11_433Asn-Lys'};

    protein_idx=logical(ismember(input_data.COMMON,genes_to_plot{1}));
    
    clear to_plot
    %also all together for YJM975
    for i=1:length(strains_to_plot)
        strain_idx_1=ismember(input_data.strain_detail,strains_to_plot{i});
        temp_idx=logical(protein_idx.*strain_idx_1);
        %sum(temp_idx)
        to_plot{i}=input_data.abundance(temp_idx);
    end


    hold on
    for j=1:length(to_plot)
        v_mean(j)=mean(to_plot{j},'omitnan');
    end
    bar(v_mean)
    for j=1:length(to_plot)
        scatter(j*ones(length(to_plot{j}),1),to_plot{j},25,'k','filled')
    end
    xticks(1:4)
    xtickangle(45)
    temp_labels={'WT',...
        '122014T-C',...
        '433Asn-Lys',...
        '122014T-C-433Asn-Lys'};
    xticklabels(temp_labels)
    title(genes_to_plot{1})
    ylim([0 6e5])
    for i=2:4
        [h p]=ttest2(to_plot{1},to_plot{i});
        plot([1 i],[5e5+2e4*i 5e5+2e4*i],'-k')
        text((1+i)/2,5e5+2e4*i+1e4,num2str(p))
    end

end


