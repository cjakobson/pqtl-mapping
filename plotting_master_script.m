%produce various plots for pQTL manuscript

clear


code_directory='/Users/cjakobson/Documents/GitHub/pqtl-mapping/';
dependency_directory='/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/pqtl-mapping-dependencies/';
output_directory='/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/manuscript-plots/';


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









