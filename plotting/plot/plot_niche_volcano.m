function []= plot_niche_volcano(niche_to_use,dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    niche_input=readtable([dependency_directory '1K_common_annotated_niche_all.csv']);

    niche_names={'Baking','Dairy','Ferm','Human','Other','Plant','Soil','Waste'};
    
    
    mapping_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv'],...
        'FileType','text');
    
    chr_idx={'I','II','III','IV','V','VI','VII','VIII','IX','X',...
        'XI','XII','XIII','XIV','XV','XVI'};
    
    for i=1:height(niche_input)
    
        temp_chr=find(ismember(chr_idx,niche_input.chr{i}));
        temp_pos=niche_input.pos(i);
    
        temp_idx=find((mapping_data.chr==temp_chr).*(mapping_data.pos==temp_pos));
    
        if ~isempty(temp_idx)
            
            niche_input.is_pqtl(i)=1;
    
            niche_input.gene{i}=mapping_data.common1{temp_idx};
    
        else
    
            niche_input.is_pqtl(i)=0;
    
        end
    
    
    end
    
    to_plot=niche_input.is_pqtl==1;
    
    v1=log10(niche_input.niche_enrichment);
    v2=niche_input.niche_q_value;
    v3=niche_input.gene;
    v4=niche_input.niche_id;
    
    v1=v1(to_plot);
    v2=v2(to_plot);
    v3=v3(to_plot);
    v4=v4(to_plot);
    
    
    %filter by which niche is enriched
    
    %subplot(2,4,5)
    hold on
    
    temp_idx=v4==niche_to_use;
    to_plot1=abs(v1(temp_idx));
    to_plot2=v2(temp_idx);
    to_plot3=v3(temp_idx);
    
    scatter(to_plot1,to_plot2,10,'k','filled')
    axis square
    xlabel('log_{10} allele enrichment')
    ylabel('log_{10}q')
    xlim([0 1.2])
    ylim([-10 25])
    plot(xlim,[2 2],':r')
    text(0.8,-3,['n sig = ' num2str(sum(to_plot2>2))])
    text(0.8,-4,['n total = ' num2str(length(to_plot2))])
    for i=1:length(to_plot3)
        if to_plot2(i)>5
            text(to_plot1(i),to_plot2(i),to_plot3{i})
        end
    end




end