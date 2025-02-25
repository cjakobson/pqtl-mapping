function [] = plot_trait_kurtosis(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;


    load([dependency_directory 'pQTLfilename.mat'])
    load([dependency_directory 'pQTLtrait.mat'])
 

    v_lambda=nan(length(filename),1);
    v_skew=nan(length(filename),1);
    v_kurtosis=nan(length(filename),1);
    for n=1:length(filename)
        
        to_get=[dependency_directory 'all-stats/' filename{n} '.mat'];
        
        if exist(to_get)
            
            load(to_get)
            
            v1=sort(-log10(pval));
            temp_epsilon=(1/length(pval));
            v2=sort(-log10(temp_epsilon:temp_epsilon:1));
    
            v_lambda(n)=median(v1)/median(v2);
            
            v_skew(n)=skewness(phenotypes);
    
            v_kurtosis(n)=kurtosis(phenotypes);
            
        end
        
    end


    histogram(v_kurtosis-3,-2:0.5:10)
    xlim([-2 10])
    xlabel('trait kurtosis-3')
    ylabel('no. of proteins')
    axis square
    text(5,100,[num2str(sum(abs(v_kurtosis-3)<2)) ' < 2 (abs.)'])



end