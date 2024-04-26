%produce various plots for pQTL manuscript

clear

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
subplot(2,3,1)
plot_reproducibility(dependency_directory,output_directory)


%C
%parents and F6 boxplots
%from pQTLplotsForManuscript.m

%Mcr1
subplot(2,6,3)
plot_parents_f6_boxplot('YKL150W',dependency_directory,output_directory)

%Gap1
subplot(2,6,4)
plot_parents_f6_boxplot('YKR039W',dependency_directory,output_directory)


%D
%heritability against abundance
%from pQTLplotsForManuscript.m
subplot(2,3,3)
plot_heritability_abundance(dependency_directory,output_directory)


%E
%effect of locus in F6 progeny [RM/YJM/RM allele/YJM allele]
%from erg11dissection.m
subplot(2,6,7)
plot_locus_effect('YKL150W','Mcr1',6952,1e5,dependency_directory,output_directory)


%F
%fraction of variance explained
%from pQTLplotsForManuscript.m
subplot(2,3,5)
plot_heritability_explained(dependency_directory,output_directory)



%G
%pQTL rarefaction
%from pQTLplotsForManuscript.m
subplot(2,3,6)
plot_rarefaction(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_1'],'-dsvg','-r0')
print([output_directory 'figure_1_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])
%H
%volcano plot (full)
%JH did plot for figures
%from proteomicsQc.m
subplot(2,3,1)
plot_main_volcano(dependency_directory,output_directory)


%I
%couple more examples

%Erg11
subplot(2,6,3)
plot_parents_f6_boxplot('YHR007C',dependency_directory,output_directory)

%Rnr4
subplot(2,6,4)
plot_parents_f6_boxplot('YGR180C',dependency_directory,output_directory)


%J
%boxplot zoom in
%JH did plot for figures
subplot(2,3,3)
plot_zoom_volcano(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_2'],'-dsvg','-r0')
print([output_directory 'figure_1_2'],'-djpeg','-r300')




%Figure S1
figure('units','normalized','outerposition',[0 0 1 1])


%A-G
%JH did plots for figures


%H
%OD correlations
%from proteomicsQc.m

%Arg4
subplot(2,3,1)
plot_od_correlations('YHR018C',dependency_directory,output_directory)

%Aco2
subplot(2,3,2)
plot_od_correlations('YJL200C',dependency_directory,output_directory)


%I
%local vs global mapping betas
%from pQTLplotsForManuscript.m
subplot(2,3,3)
plot_local_global_beta(dependency_directory,output_directory)


%J
%mapping sensitivity
%from analyzePqtlSims.m
subplot(2,3,4)
plot_sensitivity_simulations(dependency_directory,output_directory)


%K
%parent od boxplot
%from proteomicsQc.m
subplot(2,6,9)
plot_od_boxplot(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_1'],'-dsvg','-r0')
print([output_directory 'figure_S1_1'],'-djpeg','-r300')



figure('units','normalized','outerposition',[0 0 1 1])
%L
%transgression by abundance
%from proteomicsQc.m
plot_transgression(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_2'],'-dsvg','-r0')
print([output_directory 'figure_S1_2'],'-djpeg','-r300')



figure('units','normalized','outerposition',[0 0 1 1])
%M
%nPqtls as function of parental FC
%from proteomicsQc.m
subplot(2,3,1)
plot_npqtls_fc(dependency_directory,output_directory)


%N
%nPqtls as function of transgression
%from proteomicsQc.m
subplot(2,3,2)
plot_npqtls_transgression(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_3'],'-dsvg','-r0')
print([output_directory 'figure_S1_3'],'-djpeg','-r300')



close all



%Figure 2
figure('units','normalized','outerposition',[0 0 1 1])
%A
%n/a

%B
%allele effect and 1K replication
%Odc2
subplot(2,6,1)
plot_locus_effect('YOR222W','Odc2',10727,4e4,dependency_directory,output_directory)

subplot(2,6,2)
%from pQTLs1kGenomes.m
plot_locus_effect_1K('YOR222W','Odc2',10727,3e2,dependency_directory,output_directory)


%C
%Rdl1
subplot(2,6,3)
plot_locus_effect('YOR285W','Rdl1',10838,6e4,dependency_directory,output_directory)

subplot(2,6,4)
%from pQTLs1kGenomes.m
plot_locus_effect_1K('YOR285W','Rdl1',10838,1.5e3,dependency_directory,output_directory)



%D
%overall concordance
%from pQTLs1kGenomes.m
subplot(2,6,5)
plot_1K_concordance(dependency_directory,output_directory)



%E
%mRNA vs protein effects
%from assignTagSnps.m
subplot(2,3,4)
plot_beta_ase(dependency_directory,output_directory)



%F
%MS validation
%from plotValidation.m

%Ncp1
subplot(2,6,9)
plot_validation_cis('NCP1','YDJ8525 RM11 NCP1 A>T',dependency_directory,output_directory)

%Ser2
subplot(2,6,10)
plot_validation_cis('SER2','YDJ8526 RM11 SER2 G>A',dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_1'],'-dsvg','-r0')
print([output_directory 'figure_2_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])

%G
%bubble plot
%from pQTLplotsForManuscript.m
plot_bubble_plot(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_2'],'-dsvg','-r0')
print([output_directory 'figure_2_2'],'-djpeg','-r300')


%G part 2
%target counts for above bubble plot
%from pQTLplotsForManuscript.m
figure('units','normalized','outerposition',[0 0 1 1])
subplot(6,1,1)
plot_hotspot_targets(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_3'],'-dsvg','-r0')
print([output_directory 'figure_2_3'],'-djpeg','-r300')





%Figure S2
figure('units','normalized','outerposition',[0 0 1 1])


%A
%ASE reproducibility
%from pQTLplotsForManuscript.m
plot_ase_reproducibility(dependency_directory,output_directory)


%B
%pQTLs per protein
%from pQTLplotsForManuscript.m
subplot(2,3,4)
plot_npqtls_per_protein(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S2'],'-dsvg','-r0')
print([output_directory 'figure_S2'],'-djpeg','-r300')


close all



%Figure 3

%A
%n/a

figure('units','normalized','outerposition',[0 0 1 1])
pqtl_to_plot={'IRA1','IRA2','PDE2; SRA5'};
%B
%target tSNE
%from complexCovariation.m
for i=1:length(pqtl_to_plot)
    subplot(2,3,i)
    plot_pqtl_tnse(pqtl_to_plot{i},dependency_directory,output_directory)
end



%C
%cis vs trans effect size
%from pQTLplotsForManuscript.m
subplot(2,6,7)
plot_cis_trans_effect(dependency_directory,output_directory)


%D
%cis vs cumulative trans
%from pQTLplotsForManuscript.m
subplot(2,6,8)
plot_cumulative_trans_effect(dependency_directory,output_directory)


%E
%Mcr1 stair plot
%from erg11dissection.m
subplot(2,3,5)

gene_name='YKL150W';
pqtls=[6952,642,10191,10992];

plot_pqtl_stair(gene_name,pqtls,dependency_directory,output_directory)



%F
%JH plotted for fig


%G
%n/a


%H
%sign test p values
%from pQtlSignTest.m
subplot(2,3,6)

plot_pqtl_sign_test(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_3'],'-dsvg','-r0')
print([output_directory 'figure_3'],'-djpeg','-r300')






%Figure S3
figure('units','normalized','outerposition',[0 0 1 1])


%A
%pQTLs vs microarray deletions
%from compareMicroarray.m
genes={'YBR140C','YOL081W','YOR360C'};
deletions={'ira1','ira2','pde2'};
for i=1:length(genes)
    subplot(2,3,i)
    plot_pqtl_microarray(genes{i},deletions{i},dependency_directory,output_directory)
end



%B
%n/a


%C
%MAF of pQTLs vs all segregating
%from pQTLplotsForManuscript.m
subplot(2,3,4)
plot_pqtl_maf(dependency_directory,output_directory)


%D
%
%from pQtlSignTest.m
subplot(2,3,5)
plot_sign_test_by_strain(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S3'],'-dsvg','-r0')
print([output_directory 'figure_S3'],'-djpeg','-r300')



close all


%Figure 4

figure('units','normalized','outerposition',[0 0 1 1])

%A
%cis effect size
%from pQTLplotsForManuscript.m

subplot(2,6,1)
plot_cis_effect_size(dependency_directory,output_directory)


%B
%n/a


%C
%n/a


%D
%Mcr1 cis validation
subplot(2,6,2)
plot_validation_cis('MCR1','YDJ8524 RM11 MCR1 G>A',dependency_directory,output_directory)


%E
%trans effect size
%from pQTLplotsForManuscript.m
subplot(2,6,3)
plot_trans_effect_size(dependency_directory,output_directory)


%F
%IRA2 effect on Mcr1
subplot(2,6,4)
plot_locus_effect('YKL150W','Mcr1',10191,1e5,...
    dependency_directory,output_directory)


%G
%IRA2/Mcr1 reconstruction
%from plotValidation3.m
subplot(2,6,5)
plot_validation_ira2('YKL150W','YDJ8578 RM11 IRA2 G>A',...
    dependency_directory,output_directory)



%H
%BLOSUM and FoldX pQTNs vs all segregating
subplot(2,6,7)
plot_pqtn_blosum(dependency_directory,output_directory)

subplot(2,6,8)
plot_pqtn_foldx(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_4_1'],'-dsvg','-r0')
print([output_directory 'figure_4_1'],'-djpeg','-r300')



figure('units','normalized','outerposition',[0 0 1 1])

%I
%n/a



%J
%SASA and neighbors for all segregating vs all possible
%from analyse1kSecondary_v3.m
subplot(2,8,1)
plot_structure_all_possible(1,dependency_directory,output_directory)

subplot(2,8,2)
plot_structure_all_possible(2,dependency_directory,output_directory)


%K
%same for all possible vs pQTNs vs all other segr.
%from analyse1kSecondary_v3.m
subplot(2,8,3)
plot_structure_pqtn(1,dependency_directory,output_directory)

subplot(2,8,4)
plot_structure_pqtn(2,dependency_directory,output_directory)


%L
%all possible vs 1K, subdivided by secondary structure
%from analyse1kSecondary_v3.m
subplot(2,4,5)
plot_structure_1K(1,dependency_directory,output_directory)

subplot(2,4,6)
plot_structure_1K(2,dependency_directory,output_directory)


%M
%rare vs common in 1K
%from analyse1kSecondary_v3.m
subplot(2,4,7)
plot_structure_rare_common(1,dependency_directory,output_directory)

subplot(2,4,8)
plot_structure_rare_common(2,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_4_2'],'-dsvg','-r0')
print([output_directory 'figure_4_2'],'-djpeg','-r300')



