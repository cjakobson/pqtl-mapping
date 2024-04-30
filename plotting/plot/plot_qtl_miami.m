function [] = plot_qtl_miami(condition1,condition2,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %growth QTLs
    growth_input=readtable([dependency_directory 'linearNoRad.csv']);
    growth_input(growth_input.index==0,:)=[];   %geometric
    
    condition_idx1=ismember(growth_input.condition,condition1);
    v1=growth_input.index(condition_idx1);
    v2=growth_input.pVal(condition_idx1);  %top
    
    condition_idx2=ismember(growth_input.condition,condition2);
    v3=growth_input.index(condition_idx2);
    v4=-growth_input.pVal(condition_idx2);  %bottom
    
    
    hold on
    scatter(v1,v2,v2,'k','filled')
    scatter(v3,v4,-v4,grey,'filled')
    
    ylim([-200 200])
    xlim([0 12054])
    
end
    