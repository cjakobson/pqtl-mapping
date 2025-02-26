function [] = plot_overlap_heatmap(trait_order,dependency_directory,output_directory)


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

    qtl_conditions=unique(qtl_input.condition);
    qtl_conditions=qtl_conditions(trait_order);
    
    
    %overlap matrix --  within N markers
    dist_thresh=5;

    overlap_mat=nan(length(qtl_conditions));
    overlap_count=nan(length(qtl_conditions));
    condition_count=nan(length(qtl_conditions),1);
    for i=1:length(qtl_conditions)

        temp_idx1=unique(qtl_input.index(ismember(qtl_input.condition,qtl_conditions{i})));

        condition_count(i)=length(temp_idx1);
        
        for j=1:length(qtl_conditions)

            temp_idx2=unique(qtl_input.index(ismember(qtl_input.condition,qtl_conditions{j})));

            clear v_has_overlap
            for k=1:length(temp_idx1)

                v_has_overlap(k)=logical(sum(abs(temp_idx2-temp_idx1(k))<dist_thresh));

            end

            overlap_count(i,j)=sum(v_has_overlap);
            overlap_mat(i,j)=sum(v_has_overlap)/length(v_has_overlap);

        end

    end

    %add column using pQTLs
    for i=1:length(qtl_conditions)

        temp_idx1=unique(qtl_input.index(ismember(qtl_input.condition,qtl_conditions{i})));

        temp_idx2=unique(pqtl_input.index);

        clear v_has_overlap
        for k=1:length(temp_idx1)

            v_has_overlap(k)=logical(sum(abs(temp_idx2-temp_idx1(k))<dist_thresh));

        end

        overlap_count(i,length(qtl_conditions)+1)=sum(v_has_overlap);
        overlap_mat(i,length(qtl_conditions)+1)=sum(v_has_overlap)/length(v_has_overlap);

    end

    condition_count(10)
    condition_count(11)
    
    overlap_count(10,11)
    overlap_count(10,13)
    overlap_count(11,13)
    
    mean(mean(overlap_mat(:,1:12)))
    mean(overlap_mat(:,13))

    imagesc(overlap_mat,[0 1])

    m = size(get(gcf,'colormap'),1);
    %red to blue colormap
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    c = [r g b]; 
    %colormap(flipud(c))
    colormap(c)

    %colormap copper
    axis square
    colorbar
    xticks(1:length(qtl_conditions)+1)
    xtickangle(90)
    xticklabels([qtl_conditions; {'pQTLs'}])
    yticks(1:length(qtl_conditions))
    yticklabels(qtl_conditions)



    
    
end