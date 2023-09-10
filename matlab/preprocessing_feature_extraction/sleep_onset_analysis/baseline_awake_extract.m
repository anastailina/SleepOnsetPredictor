function patientObj = baseline_awake_extract(patientObj,ifcheck)

% Extracts the EEG signal corresponding to the baseline awake period if the
% subject is not discarded
%
% Author: Anastasia, Junheng Li
%
% Inputs: 
% ifcheck: this indicates whether to do checks of sleep staging before and
% after artefact correction.

%this function basically extracts the corresponding EEG signal datapoints
%for the first continious epochs (bawline awake/onset of sleep stage of
%interest)
%% Log of code

% 21/10/2021, Code started;

%% Check staging before and after artefact correction

if (1- patientObj.if_valid_stage)
    patientObj.Sbj_discard = 1;
    disp('Subject has no useful sleep onsets or baseline identified')
    return
end
if patientObj.Sbj_discard
    disp('Subject is discarded, no analysis')
    return
end

if ifcheck
    
    baseline_awake = patientObj.baseline_awake_epcs;
    baseline_awake_noart = patientObj.baseline_awake_epcs_noart;
    
    if baseline_awake(1) ~= baseline_awake_noart(1)     % check if the indices of the sleep onset differe with and without artefacts
        patientObj.if_valid_stage = 0;
        disp('Subject has inconsistent baseline awake detection results after artefact rejection, better check')
        return
    end
end



%% Take out the relevant periods
baseline_awake = patientObj.baseline_awake_epcs;
num_cont_epc = length(baseline_awake); % number of the coninious epochs that we required in hte previous sleep-staging fuction

num_eeg_ch = size(patientObj.EEG_table,1); %define the number of EEG channels 

if sum(strcmp('Baseline_Awake_Period',patientObj.EEG_table.Properties.VariableNames)) == 0  %check whether the awake period table column has already been created for this patient
    patientObj.EEG_table = addvars(patientObj.EEG_table,cell(num_eeg_ch,1),'NewVariableNames','Baseline_Awake_Period');
    
else
    warning('Baseline awake period extraction could have already been done')
end

eeg_name = patientObj.EEG_signal_name;  %extract the eeg signal name to the eeg_name
score_len = patientObj.Scoring_Epoch; %extract the durationo the scoring epoch to the score_len variable 
Fs = patientObj.EEG_sampling_freq; %extract the sampling freqeuncy to Fs
epc_data = Fs*score_len; % define the number of datapoints used

for i = 1:num_eeg_ch %iterate through all eeg channels 
    if patientObj.EEG_table.discarded(i)    % This channel will not be used if its discarded 
        continue
    end
    awake_per = [];
   
    eeg_ch = patientObj.EEG_table.(eeg_name){i};    % Current channel EEG
    if size(eeg_ch,2) > 1 %if the size of the vector containing all EEG channels is greater than 1
        eeg_ch = eeg_ch'; %transpose the eeg channel vector
    end
    
    for j = 1:num_cont_epc %iterate through the continious epochs 
        
        epc_idx_awake = baseline_awake(j); %extract the jth index of the first  basline awake periods (containing num_cont_epc number of epochs)
        epc_awake_start = (epc_idx_awake -1)*epc_data +1; %transform the epoch indices into the indices for the datapoints corresponding to such epochs 
        awake_per = [awake_per;eeg_ch(epc_awake_start:epc_awake_start+epc_data-1)]; %extract the corresponding EEG signal for the epoch int teh awake_per variable 
        
    end
    
% the algorithm automatically concatinates teh continious epochs tgether
% for the length of num_cont_epc

    patientObj.EEG_table.Baseline_Awake_Period{i} = awake_per; %upload the extracted EEG signals for the first continious awake and sleep onset periods to the EEG table of the patient strucutre
    
end

