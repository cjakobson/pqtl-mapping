function []=plot_pqtn_phenotyping(strain_to_plot,condition_to_plot,y_lim1,y_lim2,dependency_directory,output_directory)

    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    %read scanner sgaData files

    exp_name='pqtn_phenotyping';
    start_date=20230413;
    start_time=13*60+26+42/60; 

    plate_names={'a','b','c','d','e'};

    n_scanners=3;

    condition_names={'min glc',...
        'min mal',...
        'min raf',...
        'min ethanol',...
        '4NQO',...
        'doxorubicin',...
        'fluconazole',...
        'Li',...
        'Ni',...
        'rapamycin',...
        'SDS',...
        'tebuconazole',...
        'amphotericinB',...
        'caspofungin',...
        'trifluoperazine'};

    strain_names={'6635_YJM975','6649_RM11','8427_YJM975_RPA135',...
        '8429_RM11_RPA135','8524_RM11_MCR1','8525_RM11_NCP1',...
        '8526_RM11_SER2','8527_RM11_AAT2','8528_YJM975_GCS1',...
        '8529_YJM975_IRA2','8578_RM11_IRA2','RM11_ADH1',...
        'YJM975_RPI1-1','RM11_RPI1-2','YJM975_NFS1-1',...
        'YJM975_NFS1-2'};



    %get filenames
    to_read=dir([dependency_directory exp_name '/gitter']);

    %for plots
    max_growth=3000;
    max_growth2=250;

    %get sga data
    m=1;
    clear dates times time_mins scanner plates sga_mat mat_to_process ...
        temp_time growth_mat

    for i=1:length(to_read)
        
        temp_name=to_read(i).name;

        if length(temp_name)>3

            if strcmp(temp_name((end-2):end),'dat')

                to_keep{m}=temp_name;

                %grab dates and times
                temp_str=strsplit(temp_name,'_');
                dates(m)=str2num(temp_str{1});
                times{m}=temp_str{end-2};

                %convert time to minutes and normalize dates
                temp_date=dates(m)-start_date;

                temp_hours=times{m}(1:2);
                temp_mins=times{m}(3:4);
                temp_secs=times{m}(5:6);


                time_mins(m)=str2num(temp_hours)*60+str2num(temp_mins)+...
                    str2num(temp_secs)/60+24*60*temp_date-start_time;

                temp_str2=strsplit(temp_str{end},'.');

                scanner(m)=str2num(temp_str2{1}(1));
                plates{m}=temp_str2{1}(2:end);

                %get sga data
                sga_mat{m}=readtable([dependency_directory exp_name '/gitter/' temp_name]);

                mat_to_process(:,m)=sga_mat{m}.size;

                m=m+1;

            end

        end
        
    end

        
        
        
    

    m=1;
    for k=1:n_scanners

        for i=1:length(plate_names)

            scanner_idx=scanner==k;
            plate_idx=ismember(plates,plate_names{i});
            time_idx=time_mins>0;

            idx_to_use=logical(scanner_idx.*plate_idx.*time_idx);

            temp_mat=mat_to_process(:,idx_to_use);

            temp_time{m}=time_mins(idx_to_use);

            growth_mat{m}=temp_mat;

            m=m+1;

        end

    end


    %rearrange
    for i=1:length(growth_mat)

        temp_mat=growth_mat{i};

        a1idx_base=1:2:48;
        a2idx_base=2:2:48;
        b1idx_base=49:2:96;
        b2idx_base=50:2:96;

        a1idx=[];
        a2idx=[];
        b1idx=[];
        b2idx=[];

        for k=1:16

            a1idx=[a1idx 96*(k-1)+a1idx_base];
            a2idx=[a2idx 96*(k-1)+a2idx_base];
            b1idx=[b1idx 96*(k-1)+b1idx_base];
            b2idx=[b2idx 96*(k-1)+b2idx_base];

        end

        intermediate_mat=nan(size(temp_mat));

        intermediate_mat((1:384),:)=temp_mat(a1idx,:);
        intermediate_mat((384+1):(2*384),:)=temp_mat(a2idx,:);
        intermediate_mat((2*384+1):(3*384),:)=temp_mat(b1idx,:);
        intermediate_mat((3*384+1):(4*384),:)=temp_mat(b2idx,:);

        %rearrange to 96
        aidx_base=1:2:24;
        aIdx=aidx_base;

        bidx_base=2:2:24;
        bIdx=bidx_base;

        cidx_base=25:2:48;
        cIdx=cidx_base;

        didx_base=26:2:48;
        dIdx=didx_base;


        for j=1:7
            aIdx=[aIdx aidx_base+48*j];
            bIdx=[bIdx bidx_base+48*j];
            cIdx=[cIdx cidx_base+48*j];
            dIdx=[dIdx didx_base+48*j];
        end

        for j=1:4
            reorder_mat{i}(384*(j-1)+(1:96),:)=intermediate_mat(384*(j-1)+aIdx,:);
            reorder_mat{i}(384*(j-1)+((96+1):(2*96)),:)=intermediate_mat(384*(j-1)+bIdx,:);
            reorder_mat{i}(384*(j-1)+((2*96+1):(3*96)),:)=intermediate_mat(384*(j-1)+cIdx,:);
            reorder_mat{i}(384*(j-1)+((3*96+1):(4*96)),:)=intermediate_mat(384*(j-1)+dIdx,:);
        end

    end


    growthThresh=0.75;
    %norm to YJM median

    condition_idx=find(ismember(condition_names,condition_to_plot));
    
    v_median=median(reorder_mat{condition_idx},'omitnan');
    fwhm_idx=find(v_median>(max(v_median,[],'omitnan')*growthThresh));
    fwhm_idx=fwhm_idx(1);

    v_to_use=reorder_mat{condition_idx}(:,fwhm_idx);

    mutant_idx=find(ismember(strain_names,strain_to_plot));
    
    strain_idx=1;
    temp_idx=96*(strain_idx-1)+(1:96);
    to_plot{1}=v_to_use(temp_idx);

    strain_idx=2;
    temp_idx=96*(strain_idx-1)+(1:96);
    to_plot{2}=v_to_use(temp_idx);

    strain_idx=mutant_idx;
    temp_idx=96*(strain_idx-1)+(1:96);
    to_plot{3}=v_to_use(temp_idx);

    v_temp=median(to_plot{1},'omitnan');
    for k=1:length(to_plot)
        to_plot{k}=to_plot{k}./v_temp;
    end
    
    v_temp=to_plot{1};
    for k=1:length(to_plot)
        to_plot{k}=to_plot{k}./v_temp;
        to_plot{k}(isinf(to_plot{k}))=nan;
        v_mean(k)=median(to_plot{k},'omitnan');
        v_std(k)=std(to_plot{k},[],'omitnan')./sqrt(length(to_plot{k}));
    end

    hold on
    bar(v_mean)
    errorbar(1:length(v_mean),v_mean,v_std,'k.')
    %easy_box(to_plot)
    title(condition_names{condition_idx})
    xticks(1:3)
    xtickangle(45)
    xticklabels([{'YJM975','RM11'} strain_names{mutant_idx}])
    %ylim([0.6 1.2])
    ylim([y_lim1,y_lim2])
    m=1;
    for k=1:length(to_plot)
        for l=(k+1):length(to_plot)

            [h p]=ttest2(to_plot{k},to_plot{l});

            plot([k l],1.1+0.05*(m-1)*[1 1],'k')
            text((k+l)/2,1.1+0.05*(m-1),num2str(p))

            m=m+1;

        end
    end

    

    
end