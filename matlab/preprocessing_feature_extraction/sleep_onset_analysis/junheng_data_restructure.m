
%% Script to restructure Junheng's data

% Subject group
clear all
cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All'
load fteval_EWSallsbj_20Dec.mat
sbjs_ids = sbj_remained;
num_sbjs = length(sbjs_ids);


mask_30mins = epcs_to_asleep >= 61;
sbj_remained_30mins = sbj_remained(mask_30mins);


cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All/anastasia_analysis'
sbjs_to_exclude = [];

pre_sleep_duration_mine = 30 * 10;  % 30 minutes in 6-second epochs
start_recording_duration_mine = 5 * 10;  % 5 minutes in 6-second epochs
        
pre_sleep_duration_junheng = 30*20; % 30 minutes in 6-second epochs with 3 second overlap
start_recording_duration_junheng = 5*20; % 5 minutes in 6-second epochs with second overlap

pre_sleep_duration_sleep_scoring = 30*2; % 30 minutes in 30-second epochs 


for sbj_id = 1:num_sbjs
    
    if isempty(eeg_asleepper_all{sbj_id})
        sbjs_to_exclude = [sbjs_to_exclude, sbj_remained(sbj_id)];
    elseif max(size(eeg_asleepper_all{sbj_id}.ch)) < 3
        sbjs_to_exclude = [sbjs_to_exclude, sbj_remained(sbj_id)];
    else

        % Get the sleep onset epoch (in 30-second epochs)
        sleep_onset_30sec = epcs_to_asleep(sbj_id);
        
        % Convert sleep onset index from 30-second epochs to 6-second epochs (overlap 50%)
        sleep_onset_6sec = sleep_onset_30sec * 5*2;
    
       
        % Get the start points of all epochs and features for the subjed
        all_start_points = stp_ftepcs_all{sbj_id};
        old_ft_mat_sbj = ft_sbj{sbj_id}.ck;
    
        % Find the index of the sleep onset epoch
        %[~, sleep_onset_index] = min(abs(all_start_points - sleep_onset_6sec));
    
        % 1) Get the 30 minutes of EEG before sleep onset
        pre_sleep_start_index = max(1, sleep_onset_6sec - pre_sleep_duration_junheng);
        pre_sleep_end_index = sleep_onset_6sec - 1;
        pre_sleep_start_points = all_start_points(pre_sleep_start_index:2:pre_sleep_end_index);
        pre_sleep_old_fts = old_ft_mat_sbj(pre_sleep_start_index:2:pre_sleep_end_index);
    
        % 2) Get the first 5 minutes of EEG recording
        start_recording_start_index = 1;
        start_recording_end_index = min(length(all_start_points), start_recording_duration_junheng);
         % Downsampling by 2 for non-overlapping epochs
        start_recording_start_points = all_start_points(start_recording_start_index:2:start_recording_end_index);
        start_recording_old_fts = old_ft_mat_sbj(start_recording_start_index:2:start_recording_end_index);


        sbj_data = {};
        sbj_data.eeg_mat = eeg_asleepper_all{sbj_id}.ch;
        sbj_data.stp_points = stp_ftepcs_all{sbj_id};
        sbj_data.ft_mat = ft_sbj{sbj_id};
        sbj_data.ft_dscpr = ft_dscrp;
        sbj_data.Fs = Fs;
        sbj_data.epc_to_asleep = epcs_to_asleep(sbj_id);
        sbj_data.pre_sleep_start_points = pre_sleep_start_points;
        sbj_data.pre_sleep_old_fts = pre_sleep_old_fts;
        sbj_data.start_recording_start_points = start_recording_start_points;
        sbj_data.start_recording_old_fts = start_recording_old_fts;
    
        % Save the data for this subject in a separate file
        save(['subject_', num2str(sbj_remained(sbj_id)), '.mat'], 'sbj_data');
    end 
end 

sbjs_remained_final = setdiff(sbj_remained, sbjs_to_exclude);
sbjs_remained_30mins_final = setdiff(sbj_remained_30mins, sbjs_to_exclude);