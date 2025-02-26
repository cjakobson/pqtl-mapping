function [] = plot_ira2_effect_crossplot(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    validation_input=readtable([dependency_directory '20240607_proteomics-data_validation_set3.xlsx']);
    mapping_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    validation_input(ismember(validation_input.strain_id,'QC'),:)=[];
    strain_names=unique(validation_input.strain_id);
    
    protein_names=unique(validation_input.COMMON);
    protein_names(cellfun(@isempty,protein_names))=[];
    protein_names(ismember(protein_names,'NA'))=[];


    
    
    abundance_mat=nan(length(protein_names),length(strain_names));
    for i=1:length(protein_names)

        temp_table=validation_input(ismember(validation_input.COMMON,protein_names{i}),:);

        v_abundance=temp_table.abundance;

        for j=1:length(strain_names)

            v_temp_strain=v_abundance(ismember(temp_table.strain_id,strain_names{j}));

            abundance_mat(i,j)=mean(v_temp_strain,'omitnan');

        end

    end

    %plot IRA2 effects against each other and against mapping
    wt_idx=find(ismember(strain_names,'YDJ6649'));
    mut_idx=find(ismember(strain_names,'YDJ8578'));
    rm_fc=abundance_mat(:,mut_idx)./abundance_mat(:,wt_idx);
    
    wt_idx=find(ismember(strain_names,'YDJ6635'));
    mut_idx=find(ismember(strain_names,'YDJ8529'));
    yjm_fc=abundance_mat(:,mut_idx)./abundance_mat(:,wt_idx);

    %mapping betas in same order for IRA2
    locus_idx=ismember(mapping_input.index,10191);

    %match common names
    mapping_common_name=cell(height(mapping_input),1);
    for i=1:height(mapping_input)
        
        temp_str=strsplit(mapping_input.commonName{i},';');
        mapping_common_name{i}=temp_str{1};
        
    end
    
    mapping_beta=nan(length(protein_names),1);
    for i=1:length(protein_names)

        protein_idx=ismember(mapping_common_name,protein_names{i});

        temp_idx=logical(protein_idx.*locus_idx);

        if sum(temp_idx)>0

            mapping_beta(i)=mapping_input.beta(temp_idx);

        end

    end


    proteins_to_plot=~isnan(mapping_beta);

    v1=rm_fc(proteins_to_plot);
    v2=yjm_fc(proteins_to_plot);
    hold on
    scatter(log2(v1),log2(v2),10,'k','filled')
    axis square
    xlim([-1 1])
    ylim(xlim)
    %ylim([-0.6 0.6])
    ylabel('log_2 YJM IRA2*/WT')
    xlabel('log_2 RM IRA2*/WT')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')

% 
%     %label strong hits
%     temp_idx=find(abs(mapping_beta)>0.4);
%     for i=1:length(temp_idx)
%         text(log2(v1(temp_idx(i))),v2(temp_idx(i)),protein_names{temp_idx(i)})
%     end     

    
    
end