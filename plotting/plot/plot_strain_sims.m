function [] = plot_strain_sims(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    load([dependency_directory 'strainSims/pQtlFilename.mat'])
    load([dependency_directory 'strainSims/pQtlTrait.mat'])
    load([dependency_directory 'strainSims/pQtlBeta.mat'])


    groundTruthBeta=vBeta;


    qtlThresh=3;

    for n=1:length(filename)

        load([dependency_directory 'strainSims/linearPqtl/' filename{n} '.mat'])

        discoveredBeta{n}=b_fwselection;
        qtlPos{n}=bPos(pValues>qtlThresh);

        tempIsQtn=[];
        for x=1:length(candidates)

            tempIsQtn(x)=length(candidates{x})==1;

        end

        tempIsQtn=logical(tempIsQtn);

        qtnPos{n}=posToMap(tempIsQtn);

        %record sim properties to sort later
        tempStr=strsplit(filename{n},{'strains','loci','h','_','reps'});

        for x=2:length(tempStr)

            tempNum(x)=str2num(tempStr{x});

        end

        nHaps(n)=tempNum(2);
        nCausalLoci(n)=tempNum(3);
        herit(n)=tempNum(4);
        rep(n)=tempNum(5);

    end


    for i=1:length(filename)

        trueQtlPos{i}=find(vBeta{i}~=0);

        if ~isempty(trueQtlPos{i})
            for j=1:length(qtlPos{i})

                qtlDist{i}(j)=min(abs(trueQtlPos{i}-qtlPos{i}(j)));

            end

            for j=1:length(trueQtlPos{i})

                inverseQtlDist{i}(j)=min(abs(trueQtlPos{i}(j)-qtlPos{i}));

            end

            nQtl(i)=length(qtlDist{i});
            meanQtlDist(i)=mean(qtlDist{i});

            if ~isempty(qtnPos{i})
                for j=1:length(qtnPos{i})

                    qtnDist{i}(j)=min(abs(trueQtlPos{i}-qtnPos{i}(j)));

                end

                for j=1:length(trueQtlPos{i})

                    inverseQtnDist{i}(j)=min(abs(trueQtlPos{i}(j)-qtnPos{i}));

                end

                nQtn(i)=length(qtnDist{i});
                meanQtnDist(i)=mean(qtnDist{i});
            end
        end

    end



    hapNumbers=unique(nHaps);
    locusNumbers=unique(nCausalLoci);
    heritNumbers=unique(herit);

    %plot various stats
    for i=1:length(heritNumbers)

        heritIdx=herit==heritNumbers(i);

        m=1;

        for j=1:length(hapNumbers)

            hapIdx=nHaps==hapNumbers(j);

            for k=1:length(locusNumbers)

                locusIdx=nCausalLoci==locusNumbers(k);

                tempIdx=logical(heritIdx.*hapIdx.*locusIdx);%.*nRepIdx);

                %calculate stats
                for l=1:3
                    repIdx=rep==l;

                    qtlFdr{i,j,k}(l)=1-(sum(qtlDist{logical(tempIdx.*repIdx)}<5)./...
                        sum(length(qtlDist{logical(tempIdx.*repIdx)})));
                    qtlDisc{i,j,k}(l)=sum(inverseQtlDist{logical(tempIdx.*repIdx)}<5)./...
                        sum(length(inverseQtlDist{logical(tempIdx.*repIdx)}));

                    qtnFdr{i,j,k}(l)=1-(sum(qtnDist{logical(tempIdx.*repIdx)}==0)./...
                        sum(length(qtnDist{logical(tempIdx.*repIdx)})));
                    qtnDisc{i,j,k}(l)=sum(inverseQtnDist{logical(tempIdx.*repIdx)}<5)./...
                        sum(length(inverseQtnDist{logical(tempIdx.*repIdx)}));
                end


            end

        end


    end



    hapLabels={'100 segregants','250 segregants','500 segregants','1000 segregants'};
    locusLabels={'10 causal loci','50 causal loci','100 causal loci','250 causal loci'};
    %make line graphs for comparison
    i=3;
    
    hold on
    for j=1:length(hapNumbers)
        for k=1:length(locusNumbers)
            toPlot(k)=mean(qtlDisc{i,j,k});
        end
        plot(toPlot)
    end
    axis square
    legend(hapLabels)
    title(['f_{discovered QTL}' ' heritability ' num2str(heritNumbers(i)/100)])
    xtickangle(45)
    xticks(1:4)
    xticklabels(locusLabels)
    ylim([0 1])



    



    
end





