%rarefaction plots by descending abundance
function []=plot_transgression_volcano(dependency_directory,output_directory)
    
    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    [p,var_ratio,var_labels,rm_mean,yjm_mean]=...
        calculate_transgression_expectation(dependency_directory,output_directory);
    
        
    q_val=mafdr(p,'BHFDR',true);
    %q_val=p;
    
    hold on
    v1=log2(var_ratio);
    v2=-log10(q_val);
    scatter(v1,v2,10,'k','filled')
    axis square
    xlim([-4 4])
    ylim([0 200])
    xlabel('log_2(var_{F6}/var_{sim})')
    ylabel('-log_{10}q')
    text(3,50,[num2str(sum((v1>0).*(v2>2)))])
    text(-3,50,[num2str(sum((v1<0).*(v2>2)))])
    title('transgression')


    
    
end



