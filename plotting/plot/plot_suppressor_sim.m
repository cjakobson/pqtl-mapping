%plot heritability against mean abundance
function []=plot_suppressor_sim(plot_offset,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    

    load([dependency_directory 'phasedGenotype.mat'])
    hapGenotypes=phasedGenotype;
    clear phasedGenotype

    mkdir('genotypes')

    [nHaps,nPoly]=size(hapGenotypes);

    rng(0)


    %make an example trait with N loci
    %do range of suppressors
    nLoci=50;
    nSuppressors=3;


    betaMean=0.5;   %mean and SD for underlying fitness effects
    betaSD=0.1;


    hSquared=0.8;


    causalIdx=randperm(nPoly,nLoci);

    vBeta=zeros(nPoly,1);

    for m=1:length(causalIdx)

        tempLocus=causalIdx(m);

        tempBeta=betaMean+randn*betaSD*(-1)^ceil(rand*2);     %randomly assign magnitude (uniform 0 to 0.2)
        tempBeta=tempBeta*(-1)^ceil(rand*2);    %randomly assign +/-

        vBeta(tempLocus)=tempBeta;          %for RM genotypes

    end

    %hard code the modified effect
    pivotLocus=10191;
    vBeta(pivotLocus)=2*betaMean;


    suppIdx=randperm(nPoly,nSuppressors);

    for i=0:nSuppressors

        tempGenotypes=hapGenotypes;

        if i>0

            toZero=suppIdx(1:i);

            suppHaps=logical(sum(hapGenotypes(:,toZero)==1,2));

            sum(suppHaps)

            tempGenotypes(suppHaps,pivotLocus)=0;

        end


        phenotypes{i+1}=tempGenotypes*vBeta;

        phenotypes{i+1}=(phenotypes{i+1}-mean(phenotypes{i+1},'omitnan'))./...
            std(phenotypes{i+1},[],'omitnan');

        noiseMag=1-hSquared;
        noise=noiseMag*(randn(length(phenotypes{i+1}),1)).*...
            (-1).^ceil(rand(length(phenotypes{i+1}),1)*2);

        phenotypes{i+1}=phenotypes{i+1}+noise;

    end


    m=1;
    for i=1:length(phenotypes)

        subplot(2,8,plot_offset+m)
        hold on
        toPlot{1}=phenotypes{i}(hapGenotypes(:,pivotLocus)==1);
        toPlot{2}=phenotypes{i}(hapGenotypes(:,pivotLocus)==-1);
        easyBox(toPlot)
        ylim([-2 2])
        [h p]=ttest2(toPlot{1},toPlot{2});
        text(1.5,1.5,num2str(p))
        title(['nSuppressors = ' num2str(i-1)])
        xticks(1:2)
        xtickangle(45)
        xticklabels({'RM','YJM'})
        ylabel('norm. growth')

        m=m+1;
        
    end

    
end

