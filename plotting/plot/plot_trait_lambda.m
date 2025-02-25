function [] = plot_trait_lambda(dependency_directory,output_directory)


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


    histogram(v_lambda,0.5:0.05:2)
    xlim([0.8 2])
    xlabel('\lambda_{median}')
    ylabel('no. of proteins')
    axis square
    text(1.6,20,[num2str(sum(v_lambda<1.1)) ' < 1.1'])



end