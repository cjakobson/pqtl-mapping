function [true_p_sorted, expected_p_sorted] = ...
    calculate_qq(protein_idx,dependency_directory)


    load([dependency_directory 'pQTLfilename.mat'])
   
    load([dependency_directory 'linearPqtl/' filename{protein_idx} '.mat'])

    

end

