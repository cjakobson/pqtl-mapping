function [] = plot_qq_example(protein_idx,dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;


    load([dependency_directory 'pQTLfilename.mat'])
    load([dependency_directory 'pQTLtrait.mat'])
 
    
    load([dependency_directory 'all-stats/' filename{protein_idx}])

    %subplot(2,4,1)
    hold on

    v1=sort(-log10(all_p_values));
    temp_epsilon=(1/length(all_p_values));
    v2=sort(-log10(temp_epsilon:temp_epsilon:1));
    scatter(v2,v1,10,'k','filled')
    %xlim([0 -log10(temp_epsilon)])
    %ylim([0 Inf])
    xlim([0 5])
    ylim([0 50])
    plot(xlim,xlim,'-r')
    axis square
    temp_lambda=median(v1)/median(v2);
    text(3,1,num2str(temp_lambda))
    xlabel('expected')
    ylabel('actual')
    title(filename{protein_idx})




end