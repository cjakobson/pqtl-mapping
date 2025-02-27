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

    hold on
    for i=1:length(to_plot)
        histogram(to_plot{i},-1:0.05:1,'normalization','probability')
    end
    xlim([-1 1])
    ylim([0 0.12])
    legend({'all','complexes'})
    xlabel('pairwise correlation')
    ylabel('rel. freq.')
    axis square
    [p h]=ranksum(to_plot{1},to_plot{2});
    text(0.6,0.08,num2str(p))
    
end


