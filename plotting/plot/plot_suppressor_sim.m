
function []=plot_suppressor_sim(plot_offset,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    

    load([dependency_directory 'phasedGenotype.mat'])
    hap_genotypes=phasedGenotype;
    clear phasedGenotype


    [~,n_poly]=size(hap_genotypes);

    rng(0)


    %make an example trait with N loci
    %do range of suppressors
    n_loci=50;
    n_suppressors=3;


    beta_mean=0.5;   %mean and SD for underlying fitness effects
    beta_sd=0.1;


    h_squared=0.8;


    causal_idx=randperm(n_poly,n_loci);

    v_beta=zeros(n_poly,1);

    for m=1:length(causal_idx)

        temp_locus=causal_idx(m);

        temp_beta=beta_mean+randn*beta_sd*(-1)^ceil(rand*2);     %randomly assign magnitude (uniform 0 to 0.2)
        temp_beta=temp_beta*(-1)^ceil(rand*2);    %randomly assign +/-

        v_beta(temp_locus)=temp_beta;          %for RM genotypes

    end

    %hard code the modified effect
    pivot_locus=10191;
    v_beta(pivot_locus)=2*beta_mean;


    supp_idx=randperm(n_poly,n_suppressors);

    for i=0:n_suppressors

        temp_genotypes=hap_genotypes;

        if i>0

            to_zero=supp_idx(1:i);

            supp_haps=logical(sum(hap_genotypes(:,to_zero)==1,2));

            temp_genotypes(supp_haps,pivot_locus)=0;

        end


        phenotypes{i+1}=temp_genotypes*v_beta;

        phenotypes{i+1}=(phenotypes{i+1}-mean(phenotypes{i+1},'omitnan'))./...
            std(phenotypes{i+1},[],'omitnan');

        noise_mag=1-h_squared;
        noise=noise_mag*(randn(length(phenotypes{i+1}),1)).*...
            (-1).^ceil(rand(length(phenotypes{i+1}),1)*2);

        phenotypes{i+1}=phenotypes{i+1}+noise;

    end


    m=1;
    for i=1:length(phenotypes)

        subplot(2,8,plot_offset+m)
        hold on
        to_plot{1}=phenotypes{i}(hap_genotypes(:,pivot_locus)==1);
        to_plot{2}=phenotypes{i}(hap_genotypes(:,pivot_locus)==-1);
        easy_box(to_plot)
        ylim([-2 2])
        [h p]=ttest2(to_plot{1},to_plot{2});
        text(1.5,1.5,num2str(p))
        title(['nSuppressors = ' num2str(i-1)])
        xticks(1:2)
        xtickangle(45)
        xticklabels({'RM','YJM'})
        ylabel('norm. growth')

        m=m+1;
        
    end

    
end

