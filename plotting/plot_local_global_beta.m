function []=plot_local_global_beta(dependency_directory,output_directory)
    %correlation of cis effects for local v global mapping
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    allPqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    cisPqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1_cisNearby.csv']);

    %remove OD600 term
    allPqtl(allPqtl.bPos==1,:)=[];
    cisPqtl(cisPqtl.bPos==1,:)=[];

    uniqueCisProteins=unique(cisPqtl.protein);
    clear beta1 beta2
    m=1;
    for i=1:length(uniqueCisProteins)

        tempIdx1=ismember(cisPqtl.protein,uniqueCisProteins{i});
        tempQtls1=cisPqtl.index(tempIdx1);

        tempIdx2=ismember(allPqtl.protein,uniqueCisProteins{i});

        if sum(tempIdx2)>0

            tempQtls2=allPqtl.index(tempIdx2);

            intersectQtls=intersect(tempQtls1,tempQtls2);

            if ~isempty(intersectQtls)

                for k=1:length(intersectQtls)

                    beta1(m)=cisPqtl.beta(logical(tempIdx1.*...
                        ismember(cisPqtl.index,intersectQtls(k))));

                    beta2(m)=allPqtl.beta(logical(tempIdx2.*...
                        ismember(allPqtl.index,intersectQtls(k))));

                    m=m+1;

                end

            end

        end

    end

    
    hold on
    scatter(beta1,beta2,10,'k','filled')
    xlim([-1 1])
    ylim(xlim)
    xlabel('local')
    ylabel('global')
    plot(xlim,ylim,':r')
    axis square


end


