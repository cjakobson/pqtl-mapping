%plot beta vs ASE
function []=plot_validation_cis(gene_name,mutant_to_plot,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    input_data=readtable([dependency_directory 'Segregants_validation2_annotated.xlsx']);
    %mappingData=readtable('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/odCovariate/linearPqtlOd_FDR_0.1.csv');

    strain_names=unique(input_data.background_mutation);
    protein_names=unique(input_data.GENENAME);
    protein_names(cellfun(@isempty,protein_names))=[];

    %remove flc for now and regress out OD
    input_data(ismember(input_data.condition,'SM+Ura+FLC'),:)=[];

    protein_idx=logical(ismember(input_data.COMMON,gene_name));

    v_abundance=input_data.abundance(protein_idx);
    v_conc=input_data.Concentration__g__l_(protein_idx);

    %scatter(v_conc,v_abundance,'k','filled')
    temp_model=fitlm(v_conc,v_abundance);
    yI=table2array(temp_model.Coefficients(1,1));
    m=table2array(temp_model.Coefficients(2,1));
    v_temp=v_abundance-m*v_conc;
    
    strain_idx1=ismember(input_data.background_mutation,'YDJ6635 YJM975 WT');
    tempIdx=logical(protein_idx.*strain_idx1);
    to_plot{1}=input_data.abundance(tempIdx)-m*input_data.Concentration__g__l_(tempIdx);

    strain_idx2=ismember(input_data.background_mutation,'YDJ6649 RM11 WT');
    tempIdx=logical(protein_idx.*strain_idx2);
    to_plot{2}=input_data.abundance(tempIdx)-m*input_data.Concentration__g__l_(tempIdx);
    
    strain_idx3=ismember(input_data.background_mutation,mutant_to_plot);
    tempIdx=logical(protein_idx.*strain_idx3);
    to_plot{3}=input_data.abundance(tempIdx)-m*input_data.Concentration__g__l_(tempIdx);
    
    
    hold on
    for j=1:length(to_plot)
        v_mean(j)=mean(to_plot{j});
    end
    bar(v_mean)
    for j=1:length(to_plot)
        scatter(j*ones(length(to_plot{j}),1),to_plot{j},25,'k','filled')
    end
    xticks(1:3)
    xtickangle(45)
    xticklabels({'YJM975','RM11',mutant_to_plot})
    title(gene_name)
    [h p]=ttest2(to_plot{2},to_plot{3});
    plot([2 3],[1.1*mean(v_mean(2:3)) 1.1*mean(v_mean(2:3))],'-k')
    text(2.5,1.2*mean(v_mean(2:3)),num2str(p))

    
end


