function []=plot_effect_size_correlation(dependency_directory,output_directory)


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
    
    

    %correlate distances and effect sizes
    %use all pQTNs vs F6 diploid QTNs to start
    p_qtl_input=pqtl_input(pqtl_input.isQtn==1,:);
    qtl_input=qtl_input(qtl_input.isQtn==1,:);

    clear to_plot

    m=1;
    for i=1:height(p_qtl_input)

        temp_dist=abs(p_qtl_input.index(i)-qtl_input.index);

        min_dist=min(temp_dist);

        idx_to_use=find(ismember(temp_dist,min_dist));

        for j=1:length(idx_to_use)

            to_plot{1}(m)=p_qtl_input.varExp(i);
            to_plot{2}(m)=qtl_input.varExp(idx_to_use(j));

            to_plot{3}(m)=min_dist;

            %to divide by cis/trans
            to_plot{4}(m)=p_qtl_input.dist(i);

            m=m+1;

        end        

    end

    hold on
    v1=log10(to_plot{1}(to_plot{3}==0));
    v2=log10(to_plot{2}(to_plot{3}==0));
    scatter(v1,v2,'k','filled','MarkerFaceAlpha',0.25)
    xlabel('variance explained (pQTN)')
    ylabel('variance explained (exact F6 diploid QTN)')
    [r p]=corr(v1',v2');
    axis square
    % ylim([-3 0])
    % xlim([-5 0])
    xlim([-3 0])
    ylim([-4 -1])
    plot([-3 0],[-3 0],':r')
    text(-0.5,-1,num2str(r))
    text(-0.5,-1.1,num2str(p))
    
    
end

