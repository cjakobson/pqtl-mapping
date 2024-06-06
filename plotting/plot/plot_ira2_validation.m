function [] = plot_ira2_validation(plot_offset,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    pqtl_input=readtable([dependency_directory 'Segregants_validation3_annotated_v2.csv']);
    mapping_input=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    %subset to pQTNs?
    %mapping_input=mapping_input(mapping_input.isQtn==1,:);

    pqtl_input(ismember(pqtl_input.strain,'QC'),:)=[];

    strain_names=unique(pqtl_input.background_mutation);
    protein_names=unique(pqtl_input.ORF);
    protein_names(cellfun(@isempty,protein_names))=[];
    protein_names(ismember(protein_names,'NA'))=[];

    pqtl_input(ismember(pqtl_input.strain,'QC'),:)=[];

    %regress out OD
    corrected_mat=nan(length(protein_names),length(strain_names));
    for i=4:length(protein_names)

        temp_table=pqtl_input(ismember(pqtl_input.ORF,protein_names{i}),:);
        common_names{i}=temp_table.COMMON{1};

        v_abundance=temp_table.abundance;
        v_od=cellfun(@str2num,temp_table.OD600_adj);

        temp_model=fitlm(v_od,v_abundance);
        yI=table2array(temp_model.Coefficients(1,1));
        m(i)=table2array(temp_model.Coefficients(2,1));
        v_temp=v_abundance-m(i)*v_od;

        for j=1:length(strain_names)

            v_temp_strain=v_temp(ismember(temp_table.background_mutation,strain_names{j}));

            corrected_mat(i,j)=mean(v_temp_strain,'omitnan');

        end

    end


    %plot IRA2 effects against each other and against mapping
    rm_fc=corrected_mat(:,8)./corrected_mat(:,2);
    yjm_fc=corrected_mat(:,7)./corrected_mat(:,1);


    %mapping betas in same order for IRA2
    locus_idx=ismember(mapping_input.index,10191);

    mapping_beta=nan(length(protein_names),1);
    for i=1:length(protein_names)

        protein_idx=ismember(mapping_input.protein,protein_names{i});

        temp_idx=logical(protein_idx.*locus_idx);

        if sum(temp_idx)>0

            mapping_beta(i)=mapping_input.beta(temp_idx);

        end

    end




    subplot(2,4,plot_offset+1)
    v1=yjm_fc;
    v2=mapping_beta;
    hold on
    scatter(log2(v1),v2,10,'k','filled')
    axis square
    xlim([-1 1])
    ylim([-0.6 0.6])
    xlabel('log_2 YJM975 IRA2*/WT')
    ylabel('mapping \beta')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')


    %concordant/non-concordant
    n1=sum((v1>1).*(v2>0));
    n2=sum((v1<1).*(v2>0));
    n3=sum((v1>1).*(v2<0));
    n4=sum((v1<1).*(v2<0));

    temp_table=table([n1;n2],[n3;n4],...
        'VariableNames',{'mapping RM up','mapping RM down'},'RowNames',{'ms RM up','ms RM down'});
    [h,p,stats]=fishertest(temp_table);
    text(0.8,0.2,num2str(p))       

    
    
    
    v1=rm_fc;
    v2=mapping_beta;
    
    subplot(2,4,plot_offset+2)
    hold on
    scatter(log2(v1),v2,10,'k','filled')
    axis square
    xlim([-1 1])
    ylim([-0.6 0.6])
    xlabel('log_2 RM11 IRA2*/WT')
    ylabel('mapping \beta')
    plot(xlim,[0 0],':r')
    plot([0 0],ylim,':r')


    %concordant/non-concordant
    n1=sum((v1>1).*(v2>0));
    n2=sum((v1<1).*(v2>0));
    n3=sum((v1>1).*(v2<0));
    n4=sum((v1<1).*(v2<0));


    temp_table=table([n1;n2],[n3;n4],...
        'VariableNames',{'mapping RM up','mapping RM down'},'RowNames',{'ms RM up','ms RM down'});
    [h,p,stats]=fishertest(temp_table);
    text(0.8,0.2,num2str(p))       

    
    
end