%plot heritability against mean abundance
function []=plot_heritability_explained(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    [~,~,h_squared_mean,var_exp]=...
        calculate_heritability(dependency_directory,output_directory);
    
    %also compare to recent BioBank paper
    %https://www.nature.com/articles/s41586-023-06592-6
    %biobank_data=readtable([dependency_directory '41586_2023_6592_MOESM3_ESM.xlsx'],'Sheet','ST19');

    %to_plot{3}=biobank_data.TotalHeritability_TH_;
    
    %and compare to Albert mRNA
    mrna_data=readtable([dependency_directory 'elife-35471-data4-v2.xlsx']);
    
    mrna_genes=unique(mrna_data.gene);
    
    for i=1:length(mrna_genes)
        
        to_plot{2}(i)=sum(mrna_data.var_exp(ismember(mrna_data.gene,mrna_genes{i})));
        
    end
    
    hold on
    %to_plot1=var_exp./h_squared_mean;
    to_plot{1}=var_exp;
    to_plot{1}(to_plot{1}==0)=[];
    histogram(to_plot{1},0:0.05:1,'Normalization','probability')
    histogram(to_plot{2},0:0.05:1,'Normalization','probability')
    legend({'F6','Albert mRNA'})
    xlabel('variance explained')
    xlim([0 1])
    axis square
    
    to_plot{3}=var_exp./h_squared_mean;
    to_plot{3}(to_plot{3}==0)=[];
    
    median(to_plot{1},'omitnan')
    median(to_plot{2},'omitnan')
    median(to_plot{3},'omitnan')
    
    
    %also calculate median confidence intervals for mRNA and X-pQTL
    x_pqtl_data=readtable([dependency_directory '41586_2014_BFnature12904_MOESM50_ESM.xlsx']);
    
    for i=1:height(mrna_data)
        
        temp_left=strsplit(mrna_data.CI_l{i},{':','_'});
        temp_right=strsplit(mrna_data.CI_r{i},{':','_'});
        
        conf_mrna(i)=abs(str2num(temp_left{2})-str2num(temp_right{2}));
        
    end
    
    conf_x_pqtl=abs(x_pqtl_data.x2LODIntervalLeft-x_pqtl_data.x2LODIntervalRight);
    
    median(conf_mrna)
    median(conf_x_pqtl)
    
end

