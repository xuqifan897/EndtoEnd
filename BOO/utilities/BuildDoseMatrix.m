function [M,dose_data,masks]=BuildDoseMatrix(h5file, maskfile, thresh)

verbose = 2;
PTVind = 1;

[ masks ] = open_masks( maskfile, 'xyz', verbose );
PTV = masks{PTVind}.mask;
sz = size(PTV);

for i = 1:length(masks)
    masks{i}.mask = permute(masks{i}.mask,[2,1,3]);
end

dose_data = read_dose_data_permutexy( h5file, verbose, 'all', thresh);
M = dose_data.sparsemat;

dose_data = rmfield(dose_data,'sparsemat');

% dose_data.calc_metadata = dose_data_.calc_metadata;
% dose_data.beam_metadata = dose_data_.beam_metadata;
% dose_data.column_labels = dose_data_.column_labels;
% dose_data.sparsemat_info = dose_data_.sparsemat_info;
% dose_data.sparsedata = dose_data_.sparsedata;
% dose_data.read_time_sec = dose_data_.read_time_sec;

