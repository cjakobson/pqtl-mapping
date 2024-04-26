function []=plot_trans_effect_size(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    all_pqtl_data=readtable([dependency_directory  'linearPqtlOd_FDR_0.1.csv']);
    
    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];
    
    dist_idx=all_pqtl_data.dist>1e3;
    qtn_idx=all_pqtl_data.isQtn==1;
    
    variant_labels={'protein-altering','synonymous','outside ORF'};
    
    [~,v_type]=variant_types(all_pqtl_data.variantType);
    
    clear to_plot
    for i=1:length(variant_labels)

        temp_idx=logical((v_type==i)'.*qtn_idx.*dist_idx);

        to_plot{i}=all_pqtl_data.varExp(temp_idx);

    end
    
    hold on
    easy_box(to_plot)
    ylim([0 0.05])
    %ylim([0 0.15])
    xticks(1:length(variant_labels))
    xtickangle(45)
    xticklabels(variant_labels)
    ylabel('varExp')
    title('pQTLs')
    m=1;
    for i=1:length(to_plot)
        for j=(i+1):length(to_plot)
            [p h]=ranksum(to_plot{i},to_plot{j});
            text((i+j)/2,0.005+m*0.0025,num2str(p))
            m=m+1;
        end
    end
    for i=1:length(to_plot)
        text(i,0.0025,num2str(length(to_plot{i})))
    end
    
end


