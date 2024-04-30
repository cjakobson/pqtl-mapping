function []=plot_qtl_effect(locus1,locus2,input_trait,dependency_directory,output_directory)

    load([dependency_directory 'radRemapFilename.mat'])
    load([dependency_directory 'radRemapTrait.mat'])
    
    temp_trait=trait{ismember(filename,input_trait)};

    load([dependency_directory 'phasedVLCgenotype.mat'])
    
    genotypes=phasedVLCgenotype;
    clear phasedVLCgenotype

    n_strains=length(trait{1});

    [~, n_cols]=size(genotypes);
    
    %MODEL C
    %homoRM gets 1; homoYJM gets -1
    %hets all get 1
    %final: 2*nLoci columns
    [~,temp]=size(genotypes);
    n_loci=temp/4;    

    model_genotypes=zeros(n_strains,2*n_loci);
    for i=1:n_strains
        v_genotype=[genotypes(i,1:n_loci)-genotypes(i,(n_loci+1):(2*n_loci)) genotypes(i,(2*n_loci+1):(3*n_loci))+genotypes(i,(3*n_loci+1):(4*n_loci))];
        model_genotypes(i,:)=v_genotype;
    end

    clear genotypes



    %ERG11 loci to investigate
    clear genotype_idx
    

    v1=model_genotypes(:,locus1)==-1;
    v2=model_genotypes(:,locus1)==1;
    
    v3=model_genotypes(:,locus2)==-1;
    v4=model_genotypes(:,locus2)==1;

    genotype_idx{1}=find(v1.*v3);
    genotype_idx{2}=find(v1.*v4);
    genotype_idx{3}=find(v2.*v3);
    genotype_idx{4}=find(v2.*v4);
    
    
    for i=1:length(genotype_idx)
                
        v_temp{i}=temp_trait(genotype_idx{i});

        mean_to_plot(i)=mean(v_temp{i},'omitnan');
        std_to_plot(i)=std(v_temp{i},[],'omitnan')./sqrt(length(v_temp{i}));

    end

    hold on
    bar(mean_to_plot)
    errorbar(1:length(mean_to_plot),mean_to_plot,std_to_plot,'k.')
    ylim([-0.6 0.6])
    title(input_trait)
    ylabel('normalized growth')
    
    
    

end



