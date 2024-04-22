function [f_discovered,v_bins,temp_labels]=calculate_sensitivity_simulations(dependency_directory,output_directory)

    
    load([dependency_directory 'simulations/' 'pQtlFilenameSim.mat'])
    load([dependency_directory 'simulations/' 'pQtlTraitSim.mat'])
    load([dependency_directory 'simulations/' 'pQtlBetaSim.mat'])

    
    ground_truth_beta=vBeta;


    qtl_thresh=3;

    for n=1:length(filename)

        load([dependency_directory 'simulations/' 'linearPqtlSim/' filename{n} '.mat'])

        discovered_beta{n}=b_fwselection;
        qtl_pos{n}=bPos(pValues>qtl_thresh);

        temp_is_qtn=[];
        for x=1:length(candidates)

            temp_is_qtn(x)=length(candidates{x})==1;

        end

        temp_is_qtn=logical(temp_is_qtn);

        qtn_pos{n}=posToMap(temp_is_qtn);

    end


    for i=1:length(filename)

        true_qtl_pos{i}=find(vBeta{i}~=0);

        if ~isempty(true_qtl_pos{i})
            for j=1:length(qtl_pos{i})

                qtl_dist{i}(j)=min(abs(true_qtl_pos{i}-qtl_pos{i}(j)));

            end

            for j=1:length(true_qtl_pos{i})

                inverse_qtl_dist{i}(j)=min(abs(true_qtl_pos{i}(j)-qtl_pos{i}));

            end

            if ~isempty(qtn_pos{i})
                for j=1:length(qtn_pos{i})

                    qtn_dist{i}(j)=min(abs(true_qtl_pos{i}-qtn_pos{i}(j)));

                end

                for j=1:length(true_qtl_pos{i})

                    inverse_qtn_dist{i}(j)=min(abs(true_qtl_pos{i}(j)-qtn_pos{i}));

                end
            end
        end

    end

    
    
    %plot sensitivity as a function of beta
    v_dist_all=[];
    v_beta_all=[];
    for i=201:300%1:length(filename)

        v_dist_all=[v_dist_all inverse_qtl_dist{i}];

        temp_beta=vBeta{i}(vBeta{i}~=0);
        v_beta_all=[v_beta_all; temp_beta];

    end

    v_dist_all=v_dist_all';

    discovered_thresh=5;
    was_discovered=v_dist_all>discovered_thresh;
    
    v_bins=0:0.025:0.25;
    for i=1:(length(v_bins)-1)

        temp_idx=logical((v_beta_all>v_bins(i)).*(v_beta_all<=v_bins(i+1)));
        f_discovered(i)=sum(was_discovered(temp_idx),'omitnan')/sum(temp_idx,'omitnan');
        temp_labels{i}=num2str(v_bins(i+1));
    
    end

    

end


