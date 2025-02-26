function [] = plot_pqtl_targets(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    all_pqtl_data(all_pqtl_data.index==0,:)=[];   %geometric
    %all_pqtl_data=all_pqtl_data(all_pqtl_data.isQtn==1,:);
    
    pqtl_genes=unique([all_pqtl_data.gene1; all_pqtl_data.gene2]);
    pqtl_genes=pqtl_genes(2:end);
    
    
    growth_data=readtable([dependency_directory 'linearNoRad.csv']);
    growth_data(growth_data.index==0,:)=[];   %geometric
    growth_data=growth_data(growth_data.isQtn==1,:);
    
    p_thresh=3:1:15;
    clear n_targets
    for i=1:length(pqtl_genes)
        
        temp_idx=logical(ismember(all_pqtl_data.gene1,pqtl_genes{i})+...
            ismember(all_pqtl_data.gene2,pqtl_genes{i}));
        
        n_targets(i)=sum(temp_idx);
    
        temp_idx=logical(ismember(growth_data.gene1,pqtl_genes{i})+...
            ismember(growth_data.gene2,pqtl_genes{i}));
    
        n_growth(i)=length(unique(growth_data.condition(temp_idx)));
        is_growth(i)=logical(sum(temp_idx));
        
    end
    
    [v_sorted,sort_idx]=sort(n_targets,'descend');
    
    
    n_to_use=100;
    sum(v_sorted(1:n_to_use))/sum(v_sorted)
    
    sum(is_growth(sort_idx(1:n_to_use)))/n_to_use
    sum(is_growth(sort_idx((n_to_use+1):end)))/(length(is_growth)-n_to_use)
    
    
    %wrt number of targets
    clear v1 v2
    for j=1:length(p_thresh)
    
        temp_idx=all_pqtl_data.pVal>p_thresh(j);
        pqtl_genes_to_use=unique([all_pqtl_data.gene1(temp_idx);...
            all_pqtl_data.gene2(temp_idx)]);
        pqtl_genes_to_use=pqtl_genes_to_use(2:end);
    
        clear n_targets
        for i=1:length(pqtl_genes_to_use)
            
            temp_idx=logical(ismember(all_pqtl_data.gene1,pqtl_genes_to_use{i})+...
                ismember(all_pqtl_data.gene2,pqtl_genes_to_use{i}));
            
            n_targets(i)=sum(temp_idx);
        
            temp_idx=logical(ismember(growth_data.gene1,pqtl_genes_to_use{i})+...
                ismember(growth_data.gene2,pqtl_genes_to_use{i}));
        
            n_growth(i)=length(unique(growth_data.condition(temp_idx)));
            is_growth(i)=logical(sum(temp_idx));
            
        end
    
        v1(j)=mean(n_targets);
        v2(j)=std(n_targets)/sqrt(length(n_targets));
    
    end
    
    
    %subplot(2,4,5)
    hold on
    errorbar(p_thresh,v1,v2,'.k')
    axis square
    xlabel('-log_{10}p of pQTLs')
    ylabel('no. of targets')
    xlim([0 max(p_thresh)])
    ylim([0 25])
    title('pQTL genes')





end