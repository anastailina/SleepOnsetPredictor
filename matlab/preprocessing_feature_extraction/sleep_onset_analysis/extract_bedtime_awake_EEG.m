
function sbj2 = extract_bedtime_awake_EEG(sbj2, extraction_period_length)

% This function is used to extract basic information about the
% awake EEG period at bedtime such as:
%       - EEG signal data
%       - Channels of registered pre-sleep-onset EEG 
%       - Sequence of sleep stage scores for the pre-sleep-onset period
%         period
% 
% The time between sleep-onset and the start of the bedtime should be at
% least 2*extraction_period_lenght for the patient to be considered
%
% All the information about pre-sleep-onset period is summarised in the properties of
% the patient object class: sbj2.sleep_onset (with artifact
% removal) and in sbj2.sleep_onset_noart (no artifact removal)
%
% Author: Anastasia Ilina
%
%% Function inputs:
%
% sbj2:                         the patient object for processing
% extraction_period_lenght:     duration of the pre-sleep-onset period to
%                               extract (in minutes)
%
%
%% Log of code:
%
% 11/06/2023 - Created by Anastasia Ilina 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Discard checking
if sbj2.Sbj_discard == 1
    disp(['Subject ', num2str(sbj2.patient_original_id), ' discarded, not preceding'])
    return
end

if sbj2.if_valid_stage == 0 
    disp(['No valid Sleep Onset identified in Patient  ', num2str(sbj2.patient_original_id), '. They are discarded from further analysis'])
    return 
end 

if isempty(sbj2.sleep_onset)
    disp('Sleep onset has not been yet identified. Run identify_sleep_onset function first')
end

%% Start of the analysis for sleep onset without artifact (stored in sbj2.sleep_onset)

% Extract relevant information
eeg_name = sbj2.EEG_signal_name;  %extract the eeg signal name to the eeg_name
sleep_scoring = sbj2.Scoring_labels;
score_len = sbj2.Scoring_Epoch; %extract the duration of the scoring epoch to the score_len variable
score_len_in_minutes = score_len / 60;
num_epochs_to_extract = extraction_period_length / score_len_in_minutes; % calculate how many epochs to extract 
Fs = sbj2.EEG_sampling_freq; %extract the sampling freqeuncy to Fs
epc_data = Fs*score_len; % define the number of datapoints used
num_eeg_ch = size(sbj2.EEG_table,1); %define the number of EEG channels
num_sleep_onset_types = size(sbj2.sleep_onset, 1); % extract the number of different types of sleep onset definitions used 

% define an anonomized funciton handle to check whether the column name is
% already defined in a table 
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);

% initialise the table for storage of pre-sleep-onset data if it does not already exist 

if ~isTableCol(sbj2.sleep_onset,'BedtimeAwakePeriodData')
    sbj2.sleep_onset = addvars(sbj2.sleep_onset, cell(num_sleep_onset_types,1),'NewVariableNames','BedtimeAwakePeriodData');     
    for i = 1:num_sleep_onset_types
        sbj2.sleep_onset.BedtimeAwakePeriodData{i}  = table;
    end 
end


for i = 1:num_sleep_onset_types
    % Extract the sleep scoring during the pre-sleep-onset period
    bedtime_awake_table = sbj2.sleep_onset.BedtimeAwakePeriodData{i};
    awake_epochs = sbj2.base_epcs;
    epc_start_awake = awake_epochs(1);
    epcs_awake =  epc_start_awake:epc_start_awake + num_epochs_to_extract-1;
    awake_scoring = sleep_scoring(epcs_awake);
    bedtime_awake_table.TimeDuration = extraction_period_length;
    bedtime_awake_table.SleepScoring  = {awake_scoring}; 
    
    
    bedtime_awake_EEG_data = table;
    
    for ch = 1:num_eeg_ch % iterate through all EEG channels
        if sbj2.EEG_table.discarded(ch)    % This channel will not be used if its discarded
           continue
        end
        
        % get EEG data for the current EEG channel in the loop
        eeg_ch = sbj2.EEG_table.(eeg_name){ch}; 
        if size(eeg_ch,2) > 1
           eeg_ch = eeg_ch';
        end
  
        awake_period = [];
    
        for epoch = epcs_awake % iterate through the continious epochs
            epc_awake_start = (epoch - 1)*epc_data +1; % transform the epoch indices into the indices for the datapoints corresponding to such epochs
            awake_period = [awake_period;eeg_ch(epc_awake_start:epc_awake_start+epc_data-1)]; % extract the corresponding EEG signal for the epoch int teh awake_per variable
        end
        bedtime_awake_EEG_data.Channel{ch} = sbj2.EEG_table.Channel_names{ch};
        bedtime_awake_EEG_data.EEG_data{ch} = awake_period; 

   end 
   bedtime_awake_table.EEGData = {bedtime_awake_EEG_data}; 
   sbj2.sleep_onset.BedtimeAwakePeriodData{i} = bedtime_awake_table;

end 
