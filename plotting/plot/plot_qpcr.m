function []= plot_qpcr(plot_offset,dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    %plot RT-qPCR results
    
    qpcr_data=readtable([dependency_directory '20250204_qpcr.xlsx']);
    
    strain_names={'YJM','RM','RM NCP1A-177T','RM SER2G*14A'};
    
    ncp1_data=table2array(qpcr_data(1:2,17:28));
    for i=1:length(strain_names)
        ncp1_all{i}=ncp1_data(:,3*(i-1)+(1:3));
        ncp1_array{i}=reshape(mean(ncp1_data(:,3*(i-1)+(1:3))),[],1);
        ncp1_mean(i)=mean(ncp1_array{i});
    end
    
    ser2_data=table2array(qpcr_data(8:9,17:28));
    for i=1:length(strain_names)
        ser2_all{i}=ser2_data(:,3*(i-1)+(1:3));
        ser2_array{i}=reshape(mean(ser2_data(:,3*(i-1)+(1:3))),[],1);
        ser2_mean(i)=mean(reshape(ser2_array{i},[],1));
    end
    
    clear to_plot
    to_plot{1}=[2 4];
    to_plot{2}=[2 3];
    
    %SER2 (we think is mRNA-level)
    subplot(2,8,plot_offset+1)
    hold on
    bar(2.^ser2_mean(to_plot{1}))
    ylabel('2^{\DeltaC_t} SER2 vs PFK1')
    for i=1:length(to_plot{1})
        %scatter(i*ones(3),2.^ser2_array{to_plot{1}(i)},10,'k','filled')
        for j=1:3
            scatter(i*ones(2,1)-0.2+0.1*j,2.^ser2_all{to_plot{1}(i)}(:,j),20,'k','filled')
            plot(i*ones(2,1)-0.2+0.1*j,2.^ser2_all{to_plot{1}(i)}(:,j),'k')
        end
    end
    xticks(1:length(to_plot{1}))
    xtickangle(45)
    xticklabels(strain_names(to_plot{1}))
    xlim([0 3])
    ylim([0 0.16])
    [h p]=ttest2(ser2_array{to_plot{1}(1)},ser2_array{to_plot{1}(2)});
    text(1.5,0.15,num2str(p))
    
    
    
    %NCP1 (we think is post-mRNA)
    subplot(2,8,plot_offset+2)
    hold on
    bar(2.^ncp1_mean(to_plot{2}))
    ylabel('2^{\DeltaC_t} NCP1 vs PFK1')
    for i=1:length(to_plot{2})
        %scatter(i*ones(3),2.^ncp1_array{to_plot{2}(i)},10,'k','filled')
        for j=1:3
            scatter(i*ones(2)-0.2+0.1*j,2.^ncp1_all{to_plot{2}(i)}(:,j),20,'k','filled')
            plot(i*ones(2,1)-0.2+0.1*j,2.^ncp1_all{to_plot{2}(i)}(:,j),'k')
        end
    end
    xticks(1:length(to_plot{2}))
    xtickangle(45)
    xticklabels(strain_names(to_plot{2}))
    xlim([0 3])
    ylim([0 0.016])
    [h p]=ttest2(ncp1_array{to_plot{2}(1)},ncp1_array{to_plot{2}(2)});
    text(1.5,0.015,num2str(p))



end