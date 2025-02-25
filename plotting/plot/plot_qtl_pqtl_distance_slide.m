function []=plot_qtl_pqtl_distance_slide(dependency_directory,output_directory)

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


    p_thresh=3:1:15;
    clear v1 v2 v3 v4 p
    for k=1:length(p_thresh)
    
        temp_idx=logical((pqtl_input.isQtn==1).*(pqtl_input.pVal>p_thresh(k)));
        pqtn_pos_to_use=unique(pqtl_input.index(temp_idx));
    
        for i=1:length(pqtn_pos_to_use)
        
            pqtn_dist(i)=min(abs(qtn_pos-pqtn_pos_to_use(i)));
        
        end
        
        %random sets of loci of the same size
        m=1;
        for i=1:n_perms
        
            rng(i)
            v_perm=randperm(n_loci,length(pqtn_pos_to_use));
        
            for j=1:length(v_perm)
        
                not_pqtn_dist(m)=min(abs(qtn_pos-v_perm(j)));
                m=m+1;
        
            end
        
        end
    
        v1(k)=mean(pqtn_dist);
        v2(k)=mean(not_pqtn_dist);
    
        v3(k)=std(pqtn_dist)/sqrt(n_perms);
        v4(k)=std(not_pqtn_dist)/sqrt(n_perms);
    
        [h p(k)]=kstest2(pqtn_dist,not_pqtn_dist);
    
    end
    
    %subplot(2,4,2)
    %yyaxis left
    hold on
    errorbar(p_thresh,v1,v3,'.k')
    errorbar(p_thresh,v2,v4,'.r')
    xlim([0 max(p_thresh)])
    ylim([0 7])
    axis square
    xlabel('-log_{10}p of pQTNs')
    ylabel('mean distance (markers)')
    legend({'data','permutations'},'Location','southwest')
    






end


