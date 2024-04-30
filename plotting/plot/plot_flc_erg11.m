function [] = plot_flc_erg11(dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    
    
    start_date=20230426;
    start_time=16*60+13+24/60;

    plate_names={'a','b','c','d','e'};

    n_scanners=2;

    condition_names={'no drug','flc50','flc100','flc200','flc400',...
        'no drug','flc50','flc100','flc200','flc400'};

    strain_names={'YJM975 WT','122014T>C','433Asn>Lys','122014T>C+433Asn>Lys'};



    %get filenames
    to_read=dir([dependency_directory 'erg11flc/gitter']);

    %get sga data
    m=1;

    for i=1:length(to_read)
        
        temp_name=to_read(i).name;

        if length(temp_name)>3

            if strcmp(temp_name((end-2):end),'dat')

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
                sga_mat{m}=readtable([dependency_directory 'erg11flc/gitter/' temp_name]);

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

            growth_mat{m}=temp_mat;
            m=m+1;

        end

    end


    %rearrange everything to 384 first
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

        reorder_mat{i}=nan(size(temp_mat));

        reorder_mat{i}((1:384),:)=temp_mat(a1idx,:);
        reorder_mat{i}((384+1):(2*384),:)=temp_mat(a2idx,:);
        reorder_mat{i}((2*384+1):(3*384),:)=temp_mat(b1idx,:);
        reorder_mat{i}((3*384+1):(4*384),:)=temp_mat(b2idx,:);


    end

    timeToUse=60;


    for k=1:length(timeToUse)

        for i=3%1:length(reorder_mat)

            hold on

            temp_mat=reorder_mat{i};


            clear to_plot
            for j=1:length(strain_names)

                temp_idx=(384*(j-1)+1):(384*j);
                to_plot{j}=temp_mat(temp_idx,timeToUse(k));

                mean_mat(i,j)=mean(to_plot{j},'omitnan');
                sem_mat(i,j)=std(to_plot{j},[],'omitnan');

            end
            
            v_temp=to_plot{1};
            for j=1:length(to_plot)
                to_plot{j}=to_plot{j}./v_temp;
                mean_to_plot(j)=mean(to_plot{j},'omitnan');
                sem_to_plot(j)=std(to_plot{j},[],'omitnan')./sqrt(length(to_plot{j}));
            end

            %easyBox(to_plot)
            bar(mean_to_plot)
            errorbar(1:length(mean_to_plot),mean_to_plot,sem_to_plot,sem_to_plot,'.k')
            ylim([0.95 1.2])
            title([condition_names{i} ' ' num2str(timeToUse(k))])
            xticks(1:4)
            xtickangle(45)
            xticklabels(strain_names)
            for j=2:length(to_plot)
                [h p]=ttest2(to_plot{j},to_plot{1});
                text((j+1)/2,1+j/50,num2str(p))
            end

        end

    end





end



