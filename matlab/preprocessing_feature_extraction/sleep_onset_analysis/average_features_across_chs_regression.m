function averaged_samples_cell = average_features_across_chs_regression(samp_sbj2, total_num_channels)
%
% After the features for each sample were extracted (samp_sbj2.ft_vec is 
% not empty), averages features from
% the same sample across the EEG channels 
% 
%
% All the information about pre-sleep-onset period is summarised in the properties of
% the patient object class: sbj2.sleep_onset (with artifact
% removal) and in sbj2.sleep_onset_noart (no artifact removal)
%
% Author: Anastasia Ilina
%
%% Function inputs:
%
% samp_sbj2:                    cell contatining Samples of the EEG to
%                               analyse
% EEG_chs:                      list of the EEG channels in the patient 
% 
%
%% Log of code:
%
% 18/06/2023 - Created by Anastasia Ilina 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(samp_sbj2{1}.ft_vec)
   disp('Features has not been extracted yet')
   return
end 

num_samples = length(samp_sbj2);



checked_samples = [];
sbj_id = samp_sbj2{1}.Sbj_id;
averaged_samples_cell = {};
num_averaged_samples = 0; 

for sample_idx = 1:num_samples
    channel_ft_matrix = [];
    if ismember(sample_idx, checked_samples) 
        continue
    else 
        
        sample_now = samp_sbj2{sample_idx};
        label_now = sample_now.Main_label;
        num_channel = 1;
        
        ft_quality  = sample_now.ft_quality;
        num_averaged_samples = num_averaged_samples + 1;

        averaged_sample = Sample(sbj_id);
        averaged_sample.Channel_label = {sample_now.Channel_label};
        averaged_sample.Fs = sample_now.Fs;
        averaged_sample.ft_dscrp = sample_now.ft_dscrp;
        averaged_sample.Main_label = sample_now.Main_label;


        channel_ft_matrix(1,:) = sample_now.ft_vec;
        for other_sample_idx = 1:num_samples
            if other_sample_idx == sample_idx || ismember(other_sample_idx, checked_samples) 
                continue
            else   
                other_sample = samp_sbj2{other_sample_idx};
                if other_sample.Main_label== label_now
                   num_channel = num_channel + 1;
                   checked_samples = [checked_samples; other_sample_idx];
                   ft_quality = ft_quality + other_sample.ft_quality;
                   if other_sample.ft_quality ~= 0 && ~isempty(other_sample.ft_quality)
                       channel_ft_matrix(num_channel, :) = other_sample.ft_vec;
                       averaged_sample.Channel_label{num_channel} = other_sample.Channel_label;
                   end

                   if num_channel == total_num_channels
                      averaged_sample.ft_vec = mean(channel_ft_matrix, 1, 'omitnan'); 
                      averaged_sample.ft_quality = ft_quality;
                      averaged_samples_cell{num_averaged_samples} = averaged_sample;
                   end 
                   
                end 
             
            end 
        end 
    end 
end 
