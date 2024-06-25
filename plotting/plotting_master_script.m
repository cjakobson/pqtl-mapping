%produce various plots for pQTL manuscript

clear

tic

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pqtl-mapping/'];
dependency_directory=[filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/pqtl-mapping-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/manuscript-plots/'];

addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])

%Figure 1
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a

%B
%reproducibility
%from proteomicsQc.m
subplot(2,4,1)
plot_reproducibility(dependency_directory,output_directory)


%C
%volcano
%JH did plot for figures
%from proteomicsQc.m
subplot(2,4,2)
plot_main_volcano(dependency_directory,output_directory)


%D
%parents and F6 boxplots
%from pQTLplotsForManuscript.m

%Mcr1
subplot(2,8,5)
plot_parents_f6_boxplot('YKL150W',dependency_directory,output_directory)

%Gap1
subplot(2,8,6)
plot_parents_f6_boxplot('YKR039W',dependency_directory,output_directory)

%Rnr4
subplot(2,8,7)
plot_parents_f6_boxplot('YGR180C',dependency_directory,output_directory)

%Erg11
subplot(2,8,8)
plot_parents_f6_boxplot('YHR007C',dependency_directory,output_directory)





%E
%heritability against abundance
%from pQTLplotsForManuscript.m
subplot(2,4,5)
plot_heritability_abundance(dependency_directory,output_directory)


%F
%CV SGRP vs CV F6
subplot(2,4,6)
plot_cv_sgrp(dependency_directory,output_directory)



%G
%effect of locus in F6 progeny [RM/YJM/RM allele/YJM allele]
%from erg11dissection.m
subplot(2,8,13)
plot_locus_effect('YKL150W','Mcr1',6952,1e5,dependency_directory,output_directory)


%H
%n/a

%I
%Mcr1 cis validation
subplot(2,8,14)
plot_validation_cis('MCR1','YDJ8524 RM11 MCR1 G>A',dependency_directory,output_directory)



%J
%variance explained
%from pQTLplotsForManuscript.m
subplot(2,4,8)
plot_heritability_explained(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_1'],'-dsvg','-r0')
print([output_directory 'figure_1_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])


%K
%pQTL rarefaction
%from pQTLplotsForManuscript.m
subplot(2,4,1)
plot_rarefaction(dependency_directory,output_directory)





set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_2'],'-dsvg','-r0')
print([output_directory 'figure_1_2'],'-djpeg','-r300')



%Figure S1
figure('units','normalized','outerposition',[0 0 1 1])


%A-D
%JH did plots for figures


%E
%parent od boxplot
%from proteomicsQc.m
subplot(2,8,1)
plot_od_boxplot(dependency_directory,output_directory)



%F
%heritability vs C.V.
subplot(2,4,2)
plot_heritability_cv(dependency_directory,output_directory)


%G
%OD correlations
%from proteomicsQc.m

%Arg4
subplot(2,4,3)
plot_od_correlations('YHR018C',dependency_directory,output_directory)

%Aco2
subplot(2,4,4)
plot_od_correlations('YJL200C',dependency_directory,output_directory)



%H
%local vs global mapping betas
%from pQTLplotsForManuscript.m
subplot(2,4,5)
plot_local_global_beta(dependency_directory,output_directory)


%I
%mapping sensitivity
%from analyzePqtlSims.m
subplot(2,4,6)
plot_sensitivity_simulations(dependency_directory,output_directory)



%J
%boxplot zoom in
%JH did plot for figures
subplot(2,4,7)
plot_zoom_volcano(dependency_directory,output_directory)


%K
%n pQTLs vs normalized CV
subplot(2,4,8)
plot_pqtls_norm_cv(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_1'],'-dsvg','-r0')
print([output_directory 'figure_S1_1'],'-djpeg','-r300')





close all





%Figure 2
figure('units','normalized','outerposition',[0 0 1 1])
%A
%n/a

%B
%allele effect and 1K replication
%Odc2
subplot(2,8,1)
plot_locus_effect('YOR222W','Odc2',10727,4e4,dependency_directory,output_directory)

subplot(2,8,2)
%from pQTLs1kGenomes.m
plot_locus_effect_1K('YOR222W','Odc2',10727,3e2,dependency_directory,output_directory)


%C
%Rdl1
subplot(2,8,3)
plot_locus_effect('YOR285W','Rdl1',10838,6e4,dependency_directory,output_directory)

subplot(2,8,4)
%from pQTLs1kGenomes.m
plot_locus_effect_1K('YOR285W','Rdl1',10838,1.5e3,dependency_directory,output_directory)




%D
%MS validation
%from plotValidation.m

%Ncp1
subplot(2,8,5)
plot_validation_cis('NCP1','YDJ8525 RM11 NCP1 A>T',dependency_directory,output_directory)

%Ser2
subplot(2,8,6)
plot_validation_cis('SER2','YDJ8526 RM11 SER2 G>A',dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_1'],'-dsvg','-r0')
print([output_directory 'figure_2_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])

%E
%bubble plot
%from pQTLplotsForManuscript.m
plot_bubble_plot(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_2'],'-dsvg','-r0')
print([output_directory 'figure_2_2'],'-djpeg','-r300')


%E part 2
%target counts for above bubble plot
%from pQTLplotsForManuscript.m
figure('units','normalized','outerposition',[0 0 1 1])
subplot(6,1,1)
plot_hotspot_targets(dependency_directory,output_directory)


%F
%cis vs trans effect size
%from pQTLplotsForManuscript.m
subplot(2,8,9)
plot_cis_trans_effect(dependency_directory,output_directory)


%G
%cis vs cumulative trans
%from pQTLplotsForManuscript.m
subplot(2,8,10)
plot_cumulative_trans_effect(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_3'],'-dsvg','-r0')
print([output_directory 'figure_2_3'],'-djpeg','-r300')





%Figure S2
figure('units','normalized','outerposition',[0 0 1 1])



%A
%overall concordance
%from pQTLs1kGenomes.m
subplot(2,8,1)
plot_1K_concordance(dependency_directory,output_directory)



%B
%mRNA vs protein effects
%from assignTagSnps.m
subplot(2,4,2)
plot_beta_ase(dependency_directory,output_directory)



%C
%synonymous reconstructions
%GCS1
subplot(2,8,5)
plot_validation_ira2('YDL226C','YDJ8528',...
    dependency_directory,output_directory)


%AAT2
subplot(2,8,6)
plot_validation_ira2('YLR027C','YDJ8527',...
    dependency_directory,output_directory)


%D
%pQTLs per protein
%from pQTLplotsForManuscript.m
subplot(2,4,4)
plot_npqtls_per_protein(dependency_directory,output_directory)


%E
%n/a



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S2_1'],'-dsvg','-r0')
print([output_directory 'figure_S2_1'],'-djpeg','-r300')


close all



%Figure 3

figure('units','normalized','outerposition',[0 0 1 1])


%A
%n/a


%B
%Mcr1 stair plot
%from erg11dissection.m
subplot(2,4,1)

gene_name='YKL150W';
pqtls=[6952,642,10191,10992];

plot_pqtl_stair(gene_name,pqtls,dependency_directory,output_directory)




%C
%target tSNE
%from complexCovariation.m
pqtl_to_plot={'IRA1','IRA2','PDE2; SRA5'};
gene_to_highlight='YKL150W'; %Mcr1
for i=1:length(pqtl_to_plot)
    subplot(2,4,i+1)
    plot_pqtl_tnse(pqtl_to_plot{i},gene_to_highlight,dependency_directory,output_directory)
end


%D
%n/a


%E
%sign test p values
%from pQtlSignTest.m
subplot(2,4,5)

plot_pqtl_sign_test(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_3_1'],'-dsvg','-r0')
print([output_directory 'figure_3_1'],'-djpeg','-r300')






%Figure S3
figure('units','normalized','outerposition',[0 0 1 1])


%A
%pQTLs vs microarray deletions
%from compareMicroarray.m
genes={'YBR140C','YOL081W','YOR360C'};
deletions={'ira1','ira2','pde2'};
for i=1:length(genes)
    subplot(2,4,i)
    plot_pqtl_microarray(genes{i},deletions{i},dependency_directory,output_directory)
end



%B
%MAF of pQTLs vs all segregating
%from pQTLplotsForManuscript.m
subplot(2,4,4)
plot_pqtl_maf(dependency_directory,output_directory)


%C
%JH plotted for fig


%D
%
%from pQtlSignTest.m
subplot(2,4,5)
plot_sign_test_by_strain(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S3_1'],'-dsvg','-r0')
print([output_directory 'figure_S3_1'],'-djpeg','-r300')



close all





%Figure 4

figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a


%B
%cis effect size
%from pQTLplotsForManuscript.m

subplot(2,8,1)
plot_cis_effect_size(dependency_directory,output_directory)

%C
%trans effect size
%from pQTLplotsForManuscript.m
subplot(2,8,2)
plot_trans_effect_size(dependency_directory,output_directory)


%D
%IRA2 effect on Mcr1
subplot(2,8,3)
plot_locus_effect('YKL150W','Mcr1',10191,1e5,...
    dependency_directory,output_directory)


%E
%IRA2/Mcr1 reconstruction
%from plotValidation3.m
subplot(2,8,4)
plot_validation_ira2('YKL150W','YDJ8578',...
    dependency_directory,output_directory)



%F
%BLOSUM and FoldX pQTNs vs all segregating
subplot(2,8,9)
plot_pqtn_blosum(dependency_directory,output_directory)

subplot(2,8,10)
plot_pqtn_foldx(dependency_directory,output_directory)


%G
%SASA and neighbors for all segregating vs all possible
%from analyse1kSecondary_v3.m
subplot(2,8,11)
plot_structure_all_possible(1,dependency_directory,output_directory)

subplot(2,8,12)
plot_structure_all_possible(2,dependency_directory,output_directory)


%H
%same for all possible vs pQTNs vs all other segr.
%from analyse1kSecondary_v3.m
subplot(2,8,13)
plot_structure_pqtn(1,dependency_directory,output_directory)

subplot(2,8,14)
plot_structure_pqtn(2,dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_4_1'],'-dsvg','-r0')
print([output_directory 'figure_4_1'],'-djpeg','-r300')



%Figure S4
figure('units','normalized','outerposition',[0 0 1 1])



%structure property histograms
%A

%SASA
subplot(2,4,1)
plot_structure_histogram(1,5,dependency_directory,output_directory)


%neighbors
subplot(2,4,2)
plot_structure_histogram(2,1,dependency_directory,output_directory)






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S4_1'],'-dsvg','-r0')
print([output_directory 'figure_S4_1'],'-djpeg','-r300')


close all





%Figure 5
figure('units','normalized','outerposition',[0 0 1 1])
%A
%n/a


%B
%JH made fig
genes_to_plot=readtable([dependency_directory 'glycolysis_tca_genes.csv']);
name_conversion=readtable([dependency_directory 'nameConversion.csv']);

for i=1:height(name_conversion)
    
    temp_str=strsplit(name_conversion.Var1{i},';');
    name_conversion.Var1{i}=temp_str{1};
    
end

for i=1:height(genes_to_plot)
    
    genes_to_plot.commonName{i}=name_conversion.Var2{ismember(name_conversion.Var1,...
        genes_to_plot.commonName{i})};
    
end

subplot(2,4,1)
plot_correlation_heatmap(table2array(genes_to_plot),...
    dependency_directory,output_directory)


%C
%mean of parents
subplot(2,4,2)
plot_correlation_heatmap_parents(table2array(genes_to_plot),...
    dependency_directory,output_directory)


%D
%n/a


%E
%correlations amongst complex members
%from plotCorrelations.m
subplot(2,8,5)
plot_complex_correlation(dependency_directory,output_directory)


%F
%Biogrid overlap summary
subplot(2,8,6)
plot_biogrid_overlaps(dependency_directory,output_directory)



%G
%SEC61 and Sss1
subplot(2,8,7)
plot_locus_effect('YDR086C','Sss1',8221,1.8e5,dependency_directory,output_directory)


%H
%IRA2 and Bcy1
subplot(2,8,8)
plot_locus_effect('YIL033C','Bcy1',10191,1e5,dependency_directory,output_directory)


%I
%n/a


%J
%n/a


%K
%Fre1 pseudo-volcano
%from plotCorrelations.m
subplot(2,4,5)
plot_pqtl_volcano('YLR214W',dependency_directory,output_directory)



%L
%n/a


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_5_1'],'-dsvg','-r0')
print([output_directory 'figure_5_1'],'-djpeg','-r300')



%Figure S5
figure('units','normalized','outerposition',[0 0 1 1])



%A
%cross-plot various abundances
subplot(2,4,1)
gene_name1='YFR053C';   %Hxk1
gene_name2='YGL253W';   %Hxk2
plot_abundance_correlation(gene_name1,gene_name2,...
    dependency_directory,output_directory)


subplot(2,4,2)
gene_name1='YCL040W';   %Glk1
plot_abundance_correlation(gene_name1,gene_name2,...
    dependency_directory,output_directory)


%B
%correlation in correlations
%from correlaitonsAcrossExps.m
plot_corr_across_exps(2,dependency_directory,output_directory)




%C
%ATP synthase inset
%from complexCovaration.m
gene_list={'YBL099W','YJR121W','YBR039W','YPL078C','YKL016C',...
    'YPL271W','YDR377W','YML081C-A','YOL077W-A','YPR020W'};
subplot(2,4,5)
plot_correlation_heatmap(gene_list,dependency_directory,output_directory)


%C
%some other cross-plots
subplot(2,4,6)
gene_name1='YNL037C';   %Idh1
gene_name2='YOR136W';   %Idh2
plot_abundance_correlation(gene_name1,gene_name2,...
    dependency_directory,output_directory)


subplot(2,4,7)
gene_name1='YIL125W';   %Kgd1
gene_name2='YDR148C';   %Kgd2
plot_abundance_correlation(gene_name1,gene_name2,...
    dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S5_1'],'-dsvg','-r0')
print([output_directory 'figure_S5_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])



%E
%coexpression vs SWATH correlation quintiles
%from plotCorrelations.m
subplot(2,4,1)
plot_string_quintiles(4,dependency_directory,output_directory)

%F
%cellmap vs SWATH correlation quintiles
%from plotCorrelations.m
subplot(2,4,2)
plot_string_quintiles(9,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S5_2'],'-dsvg','-r0')
print([output_directory 'figure_S5_2'],'-djpeg','-r300')



close all




%Figure 6
figure('units','normalized','outerposition',[0 0 1 1])


%A
%ERG11 phenotypic mapping
%from erg11dissection.m
locus1=4972;  %Lys433Asn
locus2=4974;  %122014T>C
subplot(2,8,1)
plot_qtl_effect(locus1,locus2,'48h fluconazole_100uM-rad',dependency_directory,output_directory)



%B
%ERG11 MS validation
%from plotValidation1.m
subplot(2,8,2)
plot_validation_erg11(dependency_directory,output_directory)


%C
%ERG11 CRISPEY phenotype validation
%from readSgaDataErg11FlcGrad.m
subplot(2,8,3)
plot_flc_erg11(dependency_directory,output_directory)



%D
%NCP1 QTN vs pQTN
subplot(2,4,3)
plot_pqtn_qtn_scores_overlap(2,'YHR042W',5049,dependency_directory,output_directory)



%E
%NCP1 phenotype
subplot(2,8,7)
plot_pqtn_phenotyping('8525_RM11_NCP1','fluconazole',0.9,1.2,dependency_directory,output_directory)


%F
%n/a


%G
%pQTN and QTN scores for IRA2
%from mcr1pQtnScores.m
plot_pqtn_qtn_scores(4,'YOL081W',10191,dependency_directory,output_directory)


%H
%phenotyping on ethanol
%from readSgaDataPQtnPhen.m
subplot(2,8,13)
plot_pqtn_phenotyping('8578_RM11_IRA2','min ethanol',0.6,1.1,dependency_directory,output_directory)

subplot(2,8,14)
plot_pqtn_phenotyping('8529_YJM975_IRA2','min ethanol',0.6,1.1,dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_6_1'],'-dsvg','-r0')
print([output_directory 'figure_6_1'],'-djpeg','-r300')







%Figure S6
figure('units','normalized','outerposition',[0 0 1 1])



%A
%NCP1 phenotype on glc
subplot(2,8,1)
plot_pqtn_phenotyping('8525_RM11_NCP1','min glc',0.9,1.2,dependency_directory,output_directory)


%B
%NCP1 validation effect on Erg11 levels
subplot(2,8,2)
plot_validation_cis('ERG11','YDJ8525 RM11 NCP1 A>T',dependency_directory,output_directory)


%C
%pseudo-miami plot for pQTLs vs QTLs
subplot(2,2,2)
plot_miami('YOL081W',dependency_directory,output_directory)




%D
%IRA2 mapping vs deletion
%from compareMicroarray.m
subplot(2,4,5)
plot_pqtl_5K('YOL081W',dependency_directory,output_directory)



%E
%IRA2 reconstruction
%from plotValidation3.m
plot_ira2_validation(5,dependency_directory,output_directory)


%F
%IRA2 phenotyping in glucose
subplot(2,8,15)
plot_pqtn_phenotyping('8578_RM11_IRA2','min glc',0.6,1.1,dependency_directory,output_directory)

subplot(2,8,16)
plot_pqtn_phenotyping('8529_YJM975_IRA2','min glc',0.6,1.1,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S6_1'],'-dsvg','-r0')
print([output_directory 'figure_S6_1'],'-djpeg','-r300')


close all





figure('units','normalized','outerposition',[0 0 1 1])
%Figure 7



%A
%growth QTL rarefaction
%from pQTLvsQTL.m
subplot(2,4,1)
plot_qtl_rarefaction(dependency_directory,output_directory)

%B
%effect size of pQTLs vs QTLs
subplot(2,8,3)
plot_pqtl_qtl_effect_size(dependency_directory,output_directory)



%C
%distance from phenotypic QTNs to pQTNs/random loci
%from pQTLvsQTL.m
subplot(2,4,3)
plot_qtl_pqtl_distance(dependency_directory,output_directory)


%D
%n/a


%E
%distance from phenotypic QTNs to pQTNs/min glc QTNs
%from pQTLvsQTL.m
subplot(2,4,4)
plot_qtl_pqtl_min_glc_distance(dependency_directory,output_directory)



%F
%example miami plot for two phenotypes
subplot(2,2,3)
plot_qtl_miami('rapamycin_5uM','tebuconazole_0.6uM',dependency_directory,output_directory)


%G
%heatmap of QTN overlap
%from pQTLvsQTL.m
subplot(2,4,7)
plot_overlap_heatmap(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_7_1'],'-dsvg','-r0')
print([output_directory 'figure_7_1'],'-djpeg','-r300')






%Figure S7
figure('units','normalized','outerposition',[0 0 1 1])

%A
%pQTN vs QTN effect size
%from pQTLvsQTL.m
subplot(2,4,1)
plot_effect_size_correlation(dependency_directory,output_directory)




%B
%sensitivity from simulations
%from analyzePqtlSims.m
subplot(2,4,2)
plot_strain_sims(dependency_directory,output_directory)


%C
%heritability cross-plot
subplot(2,4,3)
plot_heritability_parents(dependency_directory,output_directory)


%D
%ASE reproducibility
%from pQTLplotsForManuscript.m
plot_ase_reproducibility(4,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S7_1'],'-dsvg','-r0')
print([output_directory 'figure_S7_1'],'-djpeg','-r300')




close all


%old stuff from 1/S1

%transgression against abundance
% subplot(2,4,5)
% plot_transgression_abundance(dependency_directory,output_directory)


%I
%fold change vs abundance
% subplot(2,4,2)
% plot_foldchange_abundance(dependency_directory,output_directory)

% 
% 
% %N
% %total variance explained vs log2 fold change
% subplot(2,4,1)
% plot_variance_foldchange(dependency_directory,output_directory)
% 
% %O
% %total variance explained vs transgression
% subplot(2,4,2)
% plot_variance_transgression(dependency_directory,output_directory)
% 
% 
% 
% 
% %Q
% %n pQTLs vs fold change
% subplot(2,4,4)
% plot_pqtls_foldchange(dependency_directory,output_directory)
% 
% 
% %R
% %n pQTLs vs transgression
% subplot(2,4,5)
% plot_pqtls_transgression(dependency_directory,output_directory)
% 


toc









