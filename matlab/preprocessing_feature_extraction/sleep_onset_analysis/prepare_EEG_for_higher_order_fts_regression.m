function sbj2 = prepare_EEG_for_higher_order_fts_regression(sbj2, EEG_ch)
% Reformat the the extracted EEG into a format suitable for extraction of
% interchannel features 
% 
%
% Author: Anastasia Ilina
%
%% Function inputs:
%
%  sbj2:                         the patient object for processing
% 
% 
%
%% Log of code:
%
% 21/06/2023 - Created by Anastasia Ilina
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check if subject discarded
if sbj2.Sbj_discard
    disp('Subject discarded, no samples returned')
    return
end

num_chs = size(sbj2.EEG_table,1); % number of channels   
sbj_id =sbj2.patient_original_id;
channel_names = sbj2.EEG_table.Channel_names;

% Check if none of the channels are discarded


sum_discarded = 0; 
for ch = 1: num_chs
    sum_discarded = sum_discarded + sbj2.EEG_table.discarded(ch);
end 

if sum_discarded > 0 
    disp('Not enough channels for higher-order interactions')
    disp(['Sbj ', num2str(sbj_id), ' is deleted from unichannel samples'])
    sbj2.Sbj_discard = 1; % discard the sbj from further analysis;
    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/unichannel'
    filename = sample_filename(sbj_id);
    delete(filename);
    return 
end
 

num_sleep_onset_types = size(sbj2.sleep_onset, 1); % extract the number of different types of sleep onset definitions used 

high_order_data_table = table;

for n = 1:num_sleep_onset_types

    sleep_onset_epochs_table = sbj2.sleep_onset.PreSleepOnset_Epoched{n};
    
    EEG_epochs_across_chs = sleep_onset_epochs_table.EEGData;
    
    miniepochs_across_channels = {};
    
    
    
    
    for ch = 1:length(EEG_epochs_across_chs)
        EEG_data_per_channel = EEG_epochs_across_chs{ch}.data;
        miniepoch_labels_per_channel = EEG_epochs_across_chs{ch}.miniepoch_labels;
        miniepoch_notnan_per_channel = EEG_epochs_across_chs{ch}.notnan;
        miniepoch_isart_per_channel = EEG_epochs_across_chs{ch}.isart;
        for sample_idx = 1:length(EEG_data_per_channel)
            if ch == 1
                miniepoch_all_channels = struct;
            else
                miniepoch_all_channels = miniepochs_across_channels{sample_idx};
            end 
           
            EEG_sample = EEG_data_per_channel{sample_idx};
            sample_label = miniepoch_labels_per_channel(sample_idx);
            ifnotnan = miniepoch_notnan_per_channel(sample_idx);
            ifisart = miniepoch_isart_per_channel(sample_idx);
            miniepoch_all_channels.sample_data(ch,:) = EEG_sample;
            
            if ch == 1
                miniepoch_all_channels.sample_label = sample_label; 
                miniepoch_all_channels.EEG_channels = EEG_ch;
            else
                if sample_label ~= miniepoch_all_channels.sample_label
                    warning ('Sample labels (time from sleep onset) varries for the same sample across different channels. RECHECK!')
                end 
            end
            miniepoch_all_channels.ifnotnan(ch) = ifnotnan;
            miniepoch_all_channels.ifisart(ch) = ifisart; 
            
            miniepochs_across_channels{sample_idx} = miniepoch_all_channels;
        end 
    end 
    
    high_order_data_table.SleepOnsetType{n} = sbj2.sleep_onset.OnsetType{n};
    high_order_data_table.EpochedData{n} = miniepochs_across_channels;

end 


sbj2.higher_order_EEG_table = high_order_data_table;