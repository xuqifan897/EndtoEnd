parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation';
target_folder = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph';
num_patients = 6;

for i = 1:num_patients
    clearvars -except parent_folder target_folder num_patients i
    contour_annealing_file = fullfile(parent_folder, ...
        ['patient', num2str(i), '_compare'], 'contour_annealing.fig');
    contour_BOO_file = fullfile(parent_folder, ...
        ['patient', num2str(i), '_compare'], 'contour_BOO.fig');
    
    output_path_annealing = fullfile(target_folder, ...
        ['patient', num2str(i), 'ContourAnnealing.png']);
    output_path_BOO = fullfile(target_folder, ...
        ['patient', num2str(i), 'ContourBOO.png']);
    openfig(contour_annealing_file);
    saveas(gcf, output_path_annealing);
    openfig(contour_BOO_file);
    saveas(gcf, output_path_BOO);
    close all
end