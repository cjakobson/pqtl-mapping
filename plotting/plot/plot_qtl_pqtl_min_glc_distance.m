function []=plot_qtl_pqtl_min_glc_distance(dependency_directory,output_directory)

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

    
    %compare min glc only QTNs vs pQTNs
    pqtn_pos=unique(pqtl_input.index(pqtl_input.isQtn==1));

    glc_idx=ismember(qtl_input.condition,'min glc_2%');
    glc_qtn_pos=unique(qtl_input.index(logical((qtl_input.isQtn==1).*glc_idx)));

    qtn_pos=unique(qtl_input.index(logical((qtl_input.isQtn==1).*~glc_idx)));

    clear v1 v2
    for i=1:length(qtn_pos)

        v_temp=abs(qtn_pos(i)-pqtn_pos);
        v1(i)=min(v_temp);

        v_temp=abs(qtn_pos(i)-glc_qtn_pos);
        v2(i)=min(v_temp);

    end

    hold on
    v_bins=0:3:100;
    histogram(v1,v_bins,'normalization','probability')
    histogram(v2,v_bins,'normalization','probability')
    title('distance from QTNs')
    legend({'to min glc pQTNs','to min glc QTNs'})
    ylabel('relative frequency')
    xlabel('distance (markers)')
    axis square
    ylim([0 0.4])
    [h p]=kstest2(v1,v2);
    text(10,0.3,num2str(p))





end


