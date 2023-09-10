function samples_sbj = sample_gen_sleepons_regression(sbj2, EEG_chs, ifcardiac, varargin)

% This function is used to generate Sample-class objects from pre-awake and post-awake 
% periods (for nocturnal awakenings) identified in Patient-class EEG data
% to prepare for feature extraction stage
%
% Input controls the columns of EEG tables to extract, as well as the
% corresponding labels of them
%
% Author: Anastasia Ilina, Junheng Li

%% Function input: 
% sbj2: the patient object for processing
% table_idx(cell of strings): names of the columns from EEG table to extract
% labels(cell of strings): labels to assign to the generated samples
% ifcardiac: flags whether the cardiac artifacts are taken into account or
% not (currently not implemented)

%% Log:
% 20/06/2023 - Adjusted to process miniepochs for sleep onset regression 
%              from sample_gen.m by Anastaia Ilina
% 13/07/2023 - Added options for processing of data with artefacts -
%               Anastasia Ilina 

%% Code

%% Varargin

if_time_limit = 0;

if isempty(varargin)
    if_no_artefact_rejection = 0;
    if_valid_stage = sbj2.if_valid_stage;

else
    for i = 1:length(varargin)
        switch varargin{i}
            case 'No Artefact Rejection'
                if_no_artefact_rejection = 1;
                if_valid_stage = sbj2.if_valid_stage_noart;
            case 'Keep clean onset'
                if_valid_stage = sbj2.if_valid_stage;
            case 'Extraction Limit (min)'
                if_time_limit = 1;
                time_limit = varargin{i+1};
                
        end
    end
end





%% Discard check 

samples_sbj = {};
% Check if subject discarded
if sbj2.Sbj_discard
    disp('Subject discarded, no samples returned')
    return
end


if if_no_artefact_rejection
    col_name ='PreSleepOnsetDataWithoutArtReject_Epoched'; 
    sleep_onset_column = 'PreSleepOnsetDataWithoutArtReject';
   
else
    col_name = 'PreSleepOnset_Epoched';
    sleep_onset_column = 'PreSleepOnsetData';
end 



num_chs = size(sbj2.EEG_table,1); % number of channels
sbj_id = sbj2.patient_original_id;

ecg_name = sbj2.ECG_signal_name;
ecg_sbj = sbj2.ECG_table.(ecg_name){1};
ecg_Fs = sbj2.ECG_sampling_freq;

eeg_Fs = sbj2.EEG_sampling_freq;

num_sample = 0;

cell_with_samples = sbj2.sleep_onset.(col_name){1,1}.EEGData;

if if_time_limit
    % Calculate the number of samples for 30 minutes (30 minutes * 60 seconds / 6 seconds per sample)
    num_samp_time_limit = time_limit * 60 / 6;
end

for ch = 1:num_chs % iterate through all the EEG channels
    if sbj2.EEG_table.discarded(ch) == 1
        continue 
    end 
    if ~isempty(cell_with_samples{ch})


        per_data = cell_with_samples{ch}.data;
        labels_now = cell_with_samples{ch}.miniepoch_labels;
        num_samp_now = size(per_data,2);

        % Calculate the starting index for the last selected minutes of the recording
        if if_time_limit
            start_idx = max(1, num_samp_now - num_samp_30_minutes + 1);
        else
            start_idx = 1;
        end 

        for i = start_idx:num_samp_now %iterate for all of the samples for that type of EEG
            samp_obj = [];       % Clear object
            num_sample = num_sample + 1;
            samp_obj = Sample(sbj_id); %create a sample object as a sample class (input subject id)
            data = per_data(:,i); %input the EEG trace for the channel n tehre
            samp_obj.data_EEG = data{1,1};
            
            samp_obj.Fs = eeg_Fs; %input the sampling frequency in there
            samp_obj.Main_label = labels_now(i); % input the label (time to sleep onset) in there
            samp_obj.Channel_label = EEG_chs{ch};
            samp_obj.sampID = i;

            if ifcardiac      % Check whether to load ECG
                samp_obj.data_ECG = ecg_sbj;
                samp_obj.Fs_ECG = ecg_Fs;
            end

            samples_sbj{num_sample} = samp_obj;
        end

    end
end

end