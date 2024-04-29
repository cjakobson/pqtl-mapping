%plot beta vs ASE
function []=plot_validation_ira2(gene_name,mutant_to_plot,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    input_data=readtable([dependency_directory 'Segregants_validation3_annotated_v2.csv']);
    %mappingData=readtable('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/odCovariate/linearPqtlOd_FDR_0.1.csv');

    
    input_data(ismember(input_data.strain,'QC'),:)=[];

    %strain_names=unique(input_data.background_mutation);
    strain_names=unique(input_data.strain);
    protein_names=unique(input_data.ORF);
    protein_names(cellfun(@isempty,protein_names))=[];
    protein_names(ismember(protein_names,'NA'))=[];

    input_data(ismember(input_data.strain,'QC'),:)=[];

    
    %regress out OD
    corrected_mat=nan(1,length(strain_names));
    
    temp_table=input_data(ismember(input_data.ORF,gene_name),:);
    common_name=temp_table.COMMON{1};

    v_abundance=temp_table.abundance;
    v_od=cellfun(@str2num,temp_table.OD600_adj);

    temp_model=fitlm(v_od,v_abundance);
    yI=table2array(temp_model.Coefficients(1,1));
    m=table2array(temp_model.Coefficients(2,1));
    v_temp=v_abundance-m*v_od;

    for j=1:length(strain_names)

        %v_tempStrain=v_temp(ismember(temp_table.background_mutation,strain_names{j}));
        v_tempStrain=v_temp(ismember(temp_table.strain,strain_names{j}));

        corrected_mat(j)=mean(v_tempStrain,'omitnan');

    end

    v_strain=[strain_names(1:2); mutant_to_plot];
    clear to_plot v_mean v_sem
    for j=1:length(v_strain)
       
        %strainIdx=ismember(temp_table.background_mutation,v_strain{j});
        strainIdx=ismember(temp_table.strain,v_strain{j});
        
        to_plot{j}=v_temp(strainIdx);
        v_mean(j)=mean(to_plot{j},'omitnan');
        v_sem(j)=std(to_plot{j},[],'omitnan')/sqrt(length(to_plot{j}));
        
    end
    
    
    hold on
    bar(v_mean)
    ylim([0 Inf])
    xticks(1:length(v_strain))
    xtickangle(45)
    xticklabels(v_strain)
    for j=1:length(to_plot)
        scatter(j*ones(length(to_plot{j}),1),to_plot{j},10,'k','filled')
    end
    title(gene_name)
    %ylim([0 125])
    ylim([0 Inf])
    [h p]=ttest2(to_plot{2},to_plot{3});
    text(2.5,120,num2str(p))
    

    
end


