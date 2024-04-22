function []=plot_npqtls_transgression(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,orf_names]=...
        parse_raw_abundance(dependency_directory,output_directory);

    all_pqtl=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);
    
    npqtls=zeros(length(orf_names),1);
    for i=1:length(orf_names)

        npqtls(i)=sum(ismember(all_pqtl.protein,orf_names{i}));

    end
    
    [transgressive_mat]=calculate_transgression(dependency_directory,output_directory);
    
    
    hold on
    v1=abs(transgressive_mat(:,4)-transgressive_mat(:,3))./...
        abs(transgressive_mat(:,2)-transgressive_mat(:,1));
    v2=npqtls;
    scatter(v1,v2,10,'k','filled')
    axis square
    set(gca,'XScale','log')
    [r p]=corr(v1,v2,'rows','complete');
    text(100,25,num2str(r))
    text(100,20,num2str(p))
    xlabel('transgression relative to parent delta')
    ylabel('N pQTLs')




end


    