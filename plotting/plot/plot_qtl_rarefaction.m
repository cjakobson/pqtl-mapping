function [] = plot_qtl_rarefaction(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    pqtl_input=readtable([dependency_directory 'linearNoRad.csv']);

    
    qtl_conditions=unique(pqtl_input.condition);
    
    
    time_idx=ismember(pqtl_input.time,'96h');

    unique_qtls=[];
    for i=1:length(qtl_conditions)

        temp_idx=logical(ismember(pqtl_input.condition,qtl_conditions{i}).*time_idx);

        unique_qtls=unique([unique_qtls; pqtl_input.index(temp_idx)]);

        v_to_plot(i)=length(unique_qtls);

    end

    hold on
    plot(v_to_plot)
    axis square
    xticks(1:length(qtl_conditions))
    xtickangle(45)
    xticklabels(qtl_conditions)
    xlim([0 13])
    ylim([0 2500])



    
end


