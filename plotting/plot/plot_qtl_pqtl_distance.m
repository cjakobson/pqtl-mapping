function []=plot_qtl_pqtl_distance(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    %load pQTL and phenotypic mapping data
    qtl_input=readtable([dependency_directory 'linearNoRad.csv']);
    
    pqtl_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    %remove OD600 term
    pqtl_input(pqtl_input.bPos==1,:)=[];

    n_loci=12054;
    n_perms=100;

    pqtn_pos=unique(pqtl_input.index(pqtl_input.isQtn==1));
    
    qtn_pos=unique(qtl_input.index(qtl_input.isQtn==1));
    for i=1:length(pqtn_pos)

        pqtn_dist(i)=min(abs(qtn_pos-pqtn_pos(i)));

    end

    %random sets of loci of the same size
    m=1;
    for i=1:n_perms

        rng(i)
        v_perm=randperm(n_loci,length(pqtn_pos));

        for j=1:length(v_perm)

            not_pqtn_dist(m)=min(abs(qtn_pos-v_perm(j)));
            m=m+1;

        end

    end

    hold on
    v_bons=0:2:50;
    histogram(pqtn_dist,v_bons,'normalization','probability')
    histogram(not_pqtn_dist,v_bons,'normalization','probability')
    title('distance from pQTNs to F_6 diploid QTNs')
    legend({'data','permutations'})
    ylabel('relative frequency')
    xlabel('distance (markers)')
    axis square
    [h p]=kstest2(pqtn_dist,not_pqtn_dist);
    text(10,0.35,num2str(p))





end


