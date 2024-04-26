function [properties_1K,properties_sim,properties_all_segregating,properties_all_other,properties_pqtn,...
    struct_mis_1K,structure_mis_sim,mis_af] = parse_structure_analysis(dependency_directory,output_directory)
    
    
    load([dependency_directory '1002dataAnn.mat'])

    temp_idx=cellfun(@isempty,gene);

    %only need gene, position, variant type, and AF for now
    gene(temp_idx)=[];
    proteinEncoded(temp_idx)=[];
    type(temp_idx)=[];
    af(temp_idx)=[];

    %minor allele
    af(af>0.5)=1-af(af>0.5);
    af(af==0)=nan;


    structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
         'bend','turn','unstr.'};

    idx_to_use=logical(ones(length(gene),1));

    load([dependency_directory 'secondary1k.mat'],'secondary')
    load([dependency_directory 'sasa1k.mat'],'sasa')
    load([dependency_directory 'neighbors1k.mat'],'neighbors')

    

    load('simulationStructureData.mat')
    %organize data to plot systematically
    property_labels={'accessible surface area','neighbors'};

    yLim1=[0 0];
    yLim2=[250 40];

    type=type(idx_to_use);
    v_mis_1K=ismember(type,'missense_variant');

    clear type

    secondary=secondary(idx_to_use);
    v_has_struct=~cellfun(@isempty,secondary);

    [~,struct_mis_1K]=structure_types(secondary(logical(v_mis_1K.*v_has_struct)));  %for categorizing

    clear secondary

    sasa=sasa(idx_to_use);
    neighbors=neighbors(idx_to_use);

    properties_1K{1,1}=sasa(logical(v_mis_1K.*v_has_struct));
    properties_1K{1,2}=neighbors(logical(v_mis_1K.*v_has_struct));


    %allele frequencies
    af=af(idx_to_use);
    mis_af=af(logical(v_mis_1K.*v_has_struct));



    structure_mis_sim=misSecondary;   %for categorizing

    properties_sim{1,1}=misSasa;
    properties_sim{1,2}=misNeighbors;

    %all segregating and trans missense pQTNs to compare
    variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

    all_pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

    %remove OD600 term
    all_pqtl_data(all_pqtl_data.bPos==1,:)=[];


    properties_all_segregating{1}=variant_info.variantSasa(ismember(variant_info.variantType,'missense_variant'));
    properties_all_segregating{2}=variant_info.variantNeighbors(ismember(variant_info.variantType,'missense_variant'));




    qtn_idx=all_pqtl_data.isQtn==1;
    mis_idx=ismember(all_pqtl_data.variantType,'missense_variant');

    pqtn_index=unique(all_pqtl_data.index(logical(qtn_idx.*mis_idx)));

    properties_pqtn{1}=variant_info.variantSasa(pqtn_index);
    properties_pqtn{2}=variant_info.variantNeighbors(pqtn_index);




    mis_idx=ismember(variant_info.variantType,'missense_variant');
    all_other_index=variant_info.index(mis_idx);
    all_other_index(ismember(all_other_index,pqtn_index))=[];

    properties_all_other{1}=variant_info.variantSasa(all_other_index);
    properties_all_other{2}=variant_info.variantNeighbors(all_other_index);


end

