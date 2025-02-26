%simulate phenotype then perform regression

function [] = linear_mixed_model_pqtl_od(protein_name,doFineMapping)


%cd(Pathname);
mkdir('linear-pqtl');


load('phasedGenotype.mat')
genotypes=phasedGenotype;

load('pQTLtrait.mat')
load('pQTLfilename.mat')
load('f6od.mat')

traitIdx=find(ismember(filename,protein_name));

filename{traitIdx}

%append harvest ods to genotype matrix (column 1)
genotypes=[f6od' genotypes];

phenotypes=trait{traitIdx};

[nStrains nCols]=size(genotypes);
nLoci=nCols;

%zero out missing growth measurements and missing genotypes
vNoSpot=isnan(phenotypes);
vNoGenotype=sum(genotypes==0,2)==nCols;

phenotypes(vNoGenotype)=0;
phenotypes(vNoSpot)=[];
genotypes(vNoSpot,:)=[];

%size(genotypes)
%size(phenotypes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assume phenotypes are already Z-scored

modelGenotypes=genotypes;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(modelGenotypes,phenotypes,'penter',10^-3,'display','on');
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);
dof = length(bPos);
pValues = -log10(stats.PVAL(bPos));
[pValues,sortIndex] = sort(pValues,'descend');
bPos = bPos(sortIndex);
toc     %this fit takes about 25min on sherlock


%now do fine mapping

%map all variants (discard those w poor pVals later)
%just consider +/-***10*** positions for now (could incorporate chr info later,
%if need be)
tic

if doFineMapping
    %neglect geometric factors (plates, edges) now and merge back later
    posToMap=bPos(pValues>3);
    
%     %remove those too close to the end (can't map)
    posToMap=posToMap(posToMap<(nLoci-10));
    posToMap=posToMap(posToMap>10);

    %calculate residuals for fine mapping
    [~,~,r] = regress(phenotypes,[ones(length(phenotypes),1),modelGenotypes(:,bPos)]);
    
    ph2=cell(length(posToMap),1);
    
    for k=1:length(posToMap)

        position1=posToMap(k);
        upper=position1+10;
        lower=position1-10;
        
        %this code from RS routine
        for i = lower:upper
            for j = lower:upper
                
                %in RS code, x is genotypes and y is phenotypes
                [ph2{k}(i-lower+1,j-lower+1)] = ...
                    fineMappingLod_multiSite_anova(i,j,position1,modelGenotypes,...
                    b_fwselection,r);

            end
        end

    end

toc
%interpret fine mapping 

candidates=cell(length(posToMap),1);

for i=1:length(ph2)
    
    [~,candidates{i}]=qtnScore(ph2{i});
    
end

vResolved=[];
for i=1:length(candidates)
    vResolved(i)=length(candidates{i})==1;
end

fracResolved=sum(vResolved)/length(vResolved);


end



%%% Calculate percentage of variance explained by each predictor in
%%% the model (only does this for regular loci, not dominant loci)
sumR = zeros(length(bPos),1);
varianceExplained = zeros(length(bPos),1);
for i = 1:length(bPos)
    newResidual = stats.yr + b_fwselection(bPos(i))*modelGenotypes(:,bPos(i));
    sumR(i) = sum(newResidual.^2) - stats.SSresid;
end
for i = 1:length(bPos)
    varianceExplained(i) = sumR(i)/sum(sumR)*dev_fwselection;
end




% Remove variables that aren't needed that would clog up HD space for when
% we save
clear genotypes;
clear stats; clear se; clear pval; clear domB;
clear inmodel; clear inmodel2; clear inmodel3; clear domSubset; clear newResidual; 
clear history; clear phasedGenotype; clear modelGenotypes;
clear secondOrderGenotype; clear trait

% Save all the variables
save(['linear-pqtl/' filename{traitIdx} '.mat']);

end



