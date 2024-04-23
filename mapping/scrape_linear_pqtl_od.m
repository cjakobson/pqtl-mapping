%scrape data from mapping output and save table of results

clear

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)
figureCounter=1;


%load info on variants
%variantInfo=readtable('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/variantInfoStructure.csv');
variantInfo=readtable('variantInfoStructure.csv');

nameInfo=readtable('nameConversion.csv');

%cis v trans
%need to quantify nearness from gene position
%consider +/- 10kb from TES and TSS as "cis"
geneInfo=readtable('TableS4.xls');
chrArray={'I','II','III','IV','V','VI','VII','VIII','IX','X',...
    'XI','XII','XIII','XIV','XV','XVI'};
for i=1:height(geneInfo)
    
    geneName{i}=geneInfo.Name{i};
    geneStart(i)=geneInfo.SGD_Start(i);
    geneEnd(i)=geneInfo.SGD_End(i);
    if sum(ismember(chrArray,geneInfo.Chrom{i}(4:end)))>0
        geneChr(i)=find(ismember(chrArray,geneInfo.Chrom{i}(4:end)));
    end

end


fdr=0.1;
load(['pval_cutoffs_from_perm_od_fdr_' num2str(fdr) '.mat'])





o=1;
for l=1:length(vCutoff)%length(filename)
    
    pValThresh=vCutoff(l);
    
    l
    
    load('pQTLfilename.mat')
    load('pQTLtrait.mat')
    load(['linear-pqtl/' filename{l} '.mat'])
    
    %filter by pVal
    bPos=bPos(pValues>pValThresh);
    varianceExplained=varianceExplained(pValues>pValThresh);
    pValues=pValues(pValues>pValThresh)
    
    isQtn=zeros(length(posToMap),1);
    hasCandidate=zeros(length(posToMap),1);
    bestCandidate=zeros(length(posToMap),1);
    for m=1:length(isQtn)
        
        vTemp=candidates{m};
        
        if length(vTemp)==1
            if vTemp==11    %QTL confirmed as QTN
                isQtn(m)=1;
                hasCandidate(m)=1;
                bestCandidate(m)=posToMap(m)+vTemp-11;
            else            %QTN shifted relative to QTL
                hasCandidate(m)=1;
                bestCandidate(m)=posToMap(m)+vTemp-11;
                isQtn(m)=1;
            end
        end
        
    end
    
    qtnPos=posToMap(logical(isQtn));
    candidatePos=posToMap(logical(hasCandidate));
    
    qtlIsQtn=ismember(bPos,qtnPos);
    qtlHasCandidate=ismember(bPos,candidatePos);
    
    for m=1:length(bPos)
        
        vCondition{o}=filename{l};
        vCommon{o}=nameInfo.Var1(ismember(nameInfo.Var2,filename{l}));
        geneDataIdx=find(ismember(geneName,filename{l}));
        if ~isempty(geneDataIdx)
            vChr(o)=geneChr(geneDataIdx);
            vStart(o)=geneStart(geneDataIdx);
            vEnd(o)=geneEnd(geneDataIdx);
        end
        
        vBpos(o)=bPos(m);
        vPval(o)=pValues(m);
        vBeta(o)=b_fwselection(bPos(m));
        vVarExp(o)=varianceExplained(m);
        
        vIndex(o)=mod(bPos(m),12054);
        
        if qtlIsQtn(m)
            vIsQtn(o)=1;
        else
            vIsQtn(o)=0;
        end
        
        if qtlHasCandidate(m)
            vFineCandidate(o)=bestCandidate(ismember(posToMap,bPos(m)));
        else
            vFineCandidate(o)=0;
        end
        
        o=o+1;
        
    end
    
end

vFineCandidate(vFineCandidate==0)=vBpos(vFineCandidate==0);

%adjust for OD600 row
vIndex(vFineCandidate>1)=vFineCandidate(vFineCandidate>1)-1;

vIndex(vFineCandidate<=1)=0;


toOutput=table(vCondition',vCommon',vChr',vStart',vEnd',vBpos',vPval',vBeta',vVarExp',vIsQtn',vFineCandidate',...
    'VariableNames',{'protein','commonName','orfChr','orfStart','orfEnd','bPos','pVal','beta','varExp',...
    'isQtn','bestCandidate'});

%make dummy info row for geometric terms
variantInfo{height(variantInfo)+1,7}={'OD600'};

vIndex(vIndex==0)=height(variantInfo);

toOutput=[toOutput variantInfo(vIndex,:)];

%need to account for SNPs inside the gene -- set to 0
for i=1:height(toOutput)
    
    temp1=min([toOutput.orfStart(i) toOutput.orfEnd(i)]);
    temp2=max([toOutput.orfStart(i) toOutput.orfEnd(i)]);
    
    if toOutput.orfChr(i)~=toOutput.chr(i)
        toOutput.dist(i)=Inf;
    elseif (temp1<toOutput.pos(i))&&...
            (temp2>toOutput.pos(i))
        toOutput.dist(i)=0;               
    else
        toOutput.dist(i)=min(abs(toOutput.orfStart(i)-toOutput.pos(i)),...
            abs(toOutput.orfEnd(i)-toOutput.pos(i)));
    end
    
end


%convert beta to percentage of mean
abundanceData=readtable('proteinMeanStdData.csv');
for i=1:height(toOutput)
    
    tempBeta=toOutput.beta(i);
    
    tempIdx=find(ismember(abundanceData.ORF,toOutput.protein{i}));
    tempMean=abundanceData.mean(tempIdx);
    tempStd=abundanceData.std(tempIdx);
    
    toOutput.percentage(i)=abs(tempBeta*tempStd)/tempMean;
    
end


%add ASE data

%also annotate with ase information
%aseData=readtable('/Users/cjakobson/Dropbox/JaroszLab/hsp90mapping/harmonizeRnaSeqAnalysis/radAseData.csv');
aseData=readtable('radAseData.csv');



for i=1:height(toOutput)
    
    if toOutput.index(i)>0

        toOutput2(i,:)=[toOutput(i,:) aseData(toOutput.index(i),[3:6 36:51])];
        
    else
        
        toOutput2(i,:)=[toOutput(i,:) array2table(zeros(1,20),'VariableNames',...
            aseData.Properties.VariableNames([3:6 36:51]))];
        
    end
    
end





writetable(toOutput2,['linear_pqtl_od_fdr_' num2str(fdr) '.csv'])





