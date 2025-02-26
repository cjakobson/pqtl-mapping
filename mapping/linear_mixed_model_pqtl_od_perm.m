%simulate phenotype then perform regression
%do N replicates and output performance

function [] = linear_mixed_model_pqtl_od_perm(protein_name,nPerms)

%cd(Pathname);
mkdir('linear-pqtl-perm');


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

for i=1:nPerms
tic
rng(i)
scramble=phenotypes(randperm(length(phenotypes)));
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(modelGenotypes,scramble,'penter',10^-3,'display','off');
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);
dof = length(bPos);
pValues{i} = -log10(stats.PVAL(bPos));
[pValues{i},sortIndex] = sort(pValues{i},'descend');
bPos = bPos(sortIndex);
toc     %this fit takes about 25min on sherlock
end




% Remove variables that aren't needed that would clog up HD space for when
% we save
clear genotypes;
clear stats; clear se; clear pval; clear domB;
clear inmodel; clear inmodel2; clear inmodel3; clear domSubset; clear newResidual; 
clear history; clear phasedGenotype; clear modelGenotypes;
clear secondOrderGenotype; clear trait

% Save all the variables
save(['linear-pqtl-perm/' filename{traitIdx} '_perm.mat']);

end



