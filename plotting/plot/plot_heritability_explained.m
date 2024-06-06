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
    biobank_data=readtable([dependency_directory '41586_2023_6592_MOESM3_ESM.xlsx'],'Sheet','ST19');

    to_plot2=biobank_data.TotalHeritability_TH_;
    
    hold on
    %to_plot1=var_exp./h_squared_mean;
    to_plot1=var_exp;
    to_plot1(to_plot1==0)=nan;
    histogram(to_plot1,0:0.05:1,'Normalization','probability')
    histogram(to_plot2,0:0.05:1,'Normalization','probability')
    legend({'F6','biobank'})
    ylabel('frequency')
    %xlabel('fraction of H^2 explained')
    xlabel('variance explained')
    xlim([0 1])
    axis square
    
    
end

