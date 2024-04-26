function []=plot_ase_reproducibility(dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;


    %ASE plot for supp
    ase_input=readtable([dependency_directory 'variantInfoExpression.csv']);
    
    for i=1:3
        
        %cast as RM higher
        v_temp1=table2array(ase_input(:,15+(3*(i-1))));
        v_temp2=table2array(ase_input(:,17+(3*(i-1))));
        
        v_temp=v_temp1;
        v_temp(ismember(v_temp2,'YJM'))=1-v_temp1(ismember(v_temp2,'YJM'));
        
        v_ase(:,i)=v_temp;
        
    end
    
    
    
    
    m=1;
    
    for i=1:3
        
        for j=(i+1):3
            
            subplot(2,3,m)
            hold on
            scatter(v_ase(:,i),v_ase(:,j),10,'k','filled')
            axis square
            xlim([0 1])
            ylim(xlim)
            xlabel(['RM rep ' num2str(i)])
            ylabel(['RM rep ' num2str(j)])
            
            [r p]=corr(v_ase(:,i),v_ase(:,j),'rows','complete');
            text(0.7,0.1,num2str(r))
            
            m=m+1;
            
        end
        
    end




end