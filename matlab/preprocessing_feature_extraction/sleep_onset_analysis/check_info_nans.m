
function [isweird, idx_to_check, nan_mask] = check_info_nans(high_order_fts_sbj2)

idx_to_check= [];
nan_mask = {};
counter = 0;
for sample_idx= 1:length(high_order_fts_sbj2)
    sample_now = high_order_fts_sbj2{sample_idx};
    ft_vec_now = sample_now.ft_vec;
    if all(isnan(ft_vec_now))

        any_nan_EEG_ch = sum(sum(isnan(sample_now.data_EEG))) < 1536;
        
        if any_nan_EEG_ch
            counter = counter +1;
            idx_to_check = [idx_to_check; sample_idx];
            nan_mask = [nan_mask; isnan(sample_now.data_EEG)];
        end 
    end 

end 

if isempty(idx_to_check)
    isweird = 0;
else 
    isweird = 1;
end 
