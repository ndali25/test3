function DataTrans = transformData(Spec, DataRaw)

DataTrans = nan(size(DataRaw));
trans     = char(Spec.Transformation);

for i = 1:length(trans)
    if strcmp(trans(i, :), 'log')
        DataTrans(:, i) = 100 * log(DataRaw(:, i));
    
    elseif strcmp(trans(i, :), 'lin')
        DataTrans(:, i) = DataRaw(:, i);
    
    else
        error('Enter valid transformation')
        
    end
end



end