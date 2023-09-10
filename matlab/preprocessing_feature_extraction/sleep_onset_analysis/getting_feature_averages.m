

new_ft_mat_pre_sleep = sbj_data_new.new_ft_mat_pre_sleep;
num_samples = length(new_ft_mat_pre_sleep);
aperiodic_mean = NaN(num_samples, 1);

for i = 1: num_samples
     ft_sample_aperiodic = new_ft_mat_pre_sleep{i}.ft_ch(:,19);
     aperiodic_mean(i) = mean(ft_sample_aperiodic, 'omitnan');
end 

plot(aperiodic_mean)
