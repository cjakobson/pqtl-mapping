function []=plot_npqtls_fc(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names]=...
        parse_raw_abundance(dependency_directory,output_directory);

    [fold_change,p_val]=calculate_parental_mean_fc(dependency_directory,output_directory);
    
    all_pqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    npqtls=zeros(length(orf_names),1);
    for i=1:length(orf_names)

        npqtls(i)=sum(ismember(all_pqtl.protein,orf_names{i}));

    end

    hold on
    v1=abs(log2(fold_change));
    v2=npqtls;
    scatter(v1,v2,10,'k','filled')
    axis square
    [r p]=corr(v1,v2,'rows','complete');
    text(4,25,num2str(r))
    text(4,20,num2str(p))
    xlabel('log_2 parent FC')
    ylabel('N pQTLs')

end


    