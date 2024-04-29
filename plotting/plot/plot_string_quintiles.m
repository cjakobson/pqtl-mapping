function [] = plot_string_quintiles(string_idx,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [string_plot_corr,string_plot_mat,string_names]=...
        parse_string_quintiles(dependency_directory,output_directory);
    
    x_vector=string_plot_corr;
    y_vector=string_plot_mat(:,string_idx);
    %redo with dot and line plots
    x_min=min(x_vector);
    x_max=max(x_vector);
    min_max_range=x_max-x_min;

    n_bins=5;
    step_size=min_max_range/n_bins;

    %make nice matrix for boxplot
    for j=1:n_bins
        bound1=x_min+(j-1)*step_size;
        bound2=x_min+j*step_size;
        temp_idx=logical((x_vector>bound1).*(x_vector<bound2));
        sliding_array{j}=y_vector(temp_idx);
    end

    dot_and_line_plot(sliding_array)
    
    xlim([0.5 n_bins+0.5])
    axis square
    [r p]=corr(x_vector,y_vector,'rows','complete');
    if string_idx<=8
        text(n_bins,100,['r = ' num2str(r)])
        text(n_bins,50,['p = ' num2str(p)])
    else
        text(n_bins,0.02,['r = ' num2str(r)])
        text(n_bins,0.015,['p = ' num2str(p)])
    end
    xlabel('SWATH abundance correlation')
    ylabel(string_names{string_idx+3})
    
    
end

