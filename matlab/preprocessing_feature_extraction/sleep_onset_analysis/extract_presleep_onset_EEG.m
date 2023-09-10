
function sbj2 = extract_presleep_onset_EEG(sbj2, extraction_period_length, varargin)

% This function is used to extract basic information about the
% pre-sleep-onset  EEG period such as:
%       - EEG singal data
%       - Channels of registered pre-sleep-onset EEG 
%       - Sequence of sleep stage scores for the pre-sleep-onset period
%         period
% 
% The time between sleep-onset and the start of the bedtime should be at
% least 2*extraction_period_lenght for the patient to be considered
%
% All the information about pre-sleep-onset period is summarised in the properties of
% the patient object class: sbj2.sleep_onset (with artifact
% removal) and in sbj2.sleep_onset_noart (no artifxact removal)
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
% 08/06/2023 - Created by Anastasia Ilina 
% 11/06/2023 - First working draft 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read varargin 

% If we are not performing artefact rejection, perform check on
% uncleaned EEG 

if_no_artefact_rejection = 0;
if_valid_stage = sbj2.if_valid_stage;

if ~isempty(varargin)
    for i = 1:length(varargin)
        switch varargin{i}
            case 'No Artefact Rejection'
                if_no_artefact_rejection = 1;
                if_valid_stage = sbj2.if_valid_stage_noart;
            case 'Keep clean onset'
                if_valid_stage = sbj2.if_valid_stage;
                
        end
    end
end



%% Discard checking
if sbj2.Sbj_discard == 1
    disp(['Subject ', num2str(sbj2.patient_original_id), ' discarded, not preceding'])
    return
end

if if_valid_stage == 0 
    disp(['No valid Sleep Onset identified in Patient  ', num2str(sbj2.patient_original_id), '. They are discarded from further analysis'])
    sbj2.Sbj_discard = 1;
    return 
end 

if isempty(sbj2.sleep_onset)
    error('Sleep onset has not been yet identified. Run identify_sleep_onset function first') 
end

if sum(sbj2.EEG_table.discarded) ~= 0 
     disp(['Subject ', num2str(sbj2.patient_original_id), ' has less than 3 channels. They are discarded, not proceding'])
     sbj2.Sbj_discard = 1; 
     return
end 

%% Reset the time extracted if the time to sleep onset is shorter than the tiem we want to extract


if extraction_period_length > sbj2.time_to_sleep
    disp(extraction_period_length)
    extraction_period_length = sbj2.time_to_sleep;
    disp(extraction_period_length)
end



%% Start of the analysis for sleep onset without artifact (stored in sbj2.sleep_onset)

% Extract relevant information
if if_no_artefact_rejection
    eeg_name = 'Signal'
else
    eeg_name = sbj2.EEG_signal_name;  %extract the eeg signal name to the eeg_name
end 

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

if if_no_artefact_rejection
    col_name ='PreSleepOnsetDataWithoutArtReject'; 
   
else
    col_name = 'PreSleepOnsetPeriodData';
end 

    if ~isTableCol(sbj2.sleep_onset, col_name)
        sbj2.sleep_onset = addvars(sbj2.sleep_onset, cell(num_sleep_onset_types,1),'NewVariableNames', col_name);     
        for i = 1:num_sleep_onset_types
            sbj2.sleep_onset.(col_name){i}  = table;
        end 
    end


for i = 1:num_sleep_onset_types
    % Extract the sleep scoring during the pre-sleep-onset period
    pre_sleep_onset_table = sbj2.sleep_onset.(col_name){i};
    sleep_onset_epochs = sbj2.sleep_onset.OnsetEpochs{i};
    epc_start_so = sleep_onset_epochs(1);
    epcs_preso =  epc_start_so-num_epochs_to_extract:epc_start_so-1;
    pre_sleep_onset_scoring = sleep_scoring(epcs_preso);
    pre_sleep_onset_table.TimeDuration = extraction_period_length;
    pre_sleep_onset_table.SleepScoring  = {pre_sleep_onset_scoring}; 
    
    varNames_EEG = {'Channel_names', 'EEG_data'};
    sleep_onset_EEG_data = table;
    pre_so_period_matrix = [];
   for ch = 1:num_eeg_ch % iterate through all EEG channels
        if sbj2.EEG_table.discarded(ch)    % This channel will not be used if its discarded
           continue
        end
        
        % get EEG data for the current EEG channel in the loop
        eeg_ch = sbj2.EEG_table.(eeg_name){ch}; 
        if size(eeg_ch,2) > 1
           eeg_ch = eeg_ch';
        end
  
        pre_so_period = [];
    
        for epoch = epcs_preso % iterate through the continious epochs
            epc_preso_start = (epoch - 1)*epc_data +1; % transform the epoch indices into the indices for the datapoints corresponding to such epochs
            pre_so_period = [pre_so_period;eeg_ch(epc_preso_start:epc_preso_start+epc_data-1)]; % extract the corresponding EEG signal for the epoch int teh awake_per variable
        end
        sleep_onset_EEG_data.Channel{ch} = sbj2.EEG_table.Channel_names{ch};
        sleep_onset_EEG_data.EEG_data{ch} = pre_so_period; 

   end 
   pre_sleep_onset_table.EEGData = {sleep_onset_EEG_data}; 
   sbj2.sleep_onset.(col_name){i} = pre_sleep_onset_table;

end 




