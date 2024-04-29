function [] = plot_complex_correlation(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [corr_string_input_data,complex_corr]=...
        parse_complex_correlation(dependency_directory,output_directory);

    clear to_plot
    to_plot{1}=corr_string_input_data.v3;
    to_plot{2}=complex_corr;

    temp_labels={'all pairs','complexes'};

    dot_and_line_plot(to_plot)
    ylim([0 0.25])
    xlim([0.5 length(to_plot)+.5])
    xticks(1:length(to_plot))
    xtickangle(45)
    xticklabels(temp_labels)
    ylabel('pairwise correlation')

    
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(1.5,0.25,num2str(p))
    
end


