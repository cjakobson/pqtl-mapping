function [string_plot_corr,string_plot_mat,string_names]=...
    parse_string_quintiles(dependency_directory,output_directory)

    corr_string_input_data=readtable([dependency_directory 'corrData_withString_withCellmap.csv']);

    string_plot_mat=table2array(corr_string_input_data(:,4:12));
    %string_plot_mat(:,9)=real(log10(string_plot_mat(:,9)));

    string_plot_corr=table2array(corr_string_input_data(:,3));

    string_plot_mat(string_plot_mat==0)=nan;

    string_names=corr_string_input_data.Properties.VariableNames;

end



