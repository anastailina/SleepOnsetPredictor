function sbj2 = pre_sleep_onset_EEG_epoching(sbj2, miniepoch_len, overlap, extraction_period_length, varargin) 

% Break up the pre-sleep-onset  EEG period into miniepochs of chosen
% duration with a pre-defined overlap 
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
% removal) and in sbj2.sleep_onset_noart (no artifact removal)
%
% Author: Anastasia Ilina
%
%% Function inputs:
%
% sbj2:                         the patient object for processing
% epoch_len:                    duration of the epochs to break down the
%                               EEG period into
% overlap:                      proportion of overlap between epochs
%
%
%% Log of code:
%
% 11/06/2023 - Created by Anastasia Ilina 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Varargin 
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
                
        end
    end
end


%% Discard checking
if sbj2.Sbj_discard == 1
    disp(['Subject ', num2str(sbj2.patient_original_id), ' discarded, not preceding'])
    return
end


num_sleep_onset_types = size(sbj2.sleep_onset,1);
num_eeg_ch = sum(sbj2.EEG_table.discarded == 0); %define the number of EEG channels 
channels_list = [1,2,3];
channels_list = channels_list(sbj2.EEG_table.discarded == 0);

%% Reset the time extracted if the time to sleep onset is shorter than the tiem we want to extract


if extraction_period_length > sbj2.time_to_sleep

    extraction_period_length = sbj2.time_to_sleep;

end

% define an anonomized funciton handle to check whether the column name is
% already defined in a table 
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);

% initialise the table for storage of pre-sleep-onset data if it does not already exist 

if if_no_artefact_rejection
    col_name ='PreSleepOnsetDataWithoutArtReject_Epoched'; 
    sleep_onset_column = 'PreSleepOnsetDataWithoutArtReject';
   
else
    col_name = 'PreSleepOnset_Epoched';
    sleep_onset_column = 'PreSleepOnsetData';
end 


if ~isTableCol(sbj2.sleep_onset,col_name)
    sbj2.sleep_onset = addvars(sbj2.sleep_onset,cell(num_sleep_onset_types,1),'NewVariableNames',col_name);
end 

Fs = sbj2.EEG_sampling_freq; %extract the sampling freqeuncy to Fs
score_len = sbj2.Scoring_Epoch; %extract the durationo the scoring epoch to the score_len variable 

eeg_name = sbj2.EEG_signal_name;  %extract the eeg signal name to the eeg_name 
miniepoch_data = Fs*miniepoch_len; % define the number of datapoints used
jump = floor((1-overlap)*miniepoch_data); %this is by how much we 'move' when with each iteration when dividing the eeg signal into even smaller epochs 


for k = 1:num_sleep_onset_types
    

    preso_table = sbj2.sleep_onset.(sleep_onset_column){k};

    epoched_table = table;
        
    epoched_table.TimeDuration  = preso_table.TimeDuration;
    epoched_table.SleepScoring = preso_table.SleepScoring;
    all_channels = preso_table.EEGData{1,1};
    miniepochs = {}; 
    for ch_idx= 1:num_eeg_ch

        ch = channels_list(ch_idx);

        
        channel_EEG_data = all_channels.EEG_data(ch);
        channel_EEG_data = channel_EEG_data{1,1};
   

        if ~isempty(channel_EEG_data) 
            
              
            if size(channel_EEG_data ,1) > 1
               channel_EEG_data  = channel_EEG_data';
            end
              

            if overlap == 0
                len_data = size(channel_EEG_data,2);
            else
                len_data = size(channel_EEG_data,2)  - miniepoch_data;
            end 
            
            eeg_sig_ind2start = (miniepoch_data - jump) +1;
            %num_epcs_max = floor((len_data-miniepoch_data)/jump)+1; % that's
            %how junheng did it
            num_epcs_max = floor((len_data)/jump);
            distance2sleep_onset = zeros([1, num_epcs_max]);
            isart = zeros(1,num_epcs_max);    % Logic indicater of artefacts
            epcs_all_cell = cell(1, num_epcs_max);
            per_data = channel_EEG_data(eeg_sig_ind2start:end);
            epcs_all = zeros(num_epcs_max,miniepoch_data);% a matrix consisting of eachsmall epoch in each row 
            for i = 1:num_epcs_max
                distance2sleep_onset(1,i) = extraction_period_length - ((i-1)*miniepoch_len/60); % for each sample, calculate the distance to sleep onset in minutes 
                stp = 1+(i-1)*jump;
                epcs_all(i,:) = per_data(stp:stp+miniepoch_data-1); %in each row of the matri we're uploaing the datapoints from the smaller epoch
            end 
            
            if if_no_artefact_rejection
                isart = [];
            else
                % calculating the artifact threshold
                epc_rms = rms(epcs_all,2);    % Compute the root mean square for each epoch
                amp_th = 3*median(epc_rms, 'omitnan');      % Amplitude threshold defined as 3 times the mean RMS
                %artifact detection 
                for i = 1:num_epcs_max
                    if max(abs(epcs_all(i,:))) > amp_th
                        isart(i) = 1;
                        epcs_all(i,:) = NaN;
                    end
                end
            end

            miniepoch_notnan =[];
            %count NaNs
            for i = 1:num_epcs_max
                miniepoch_notnan(i) = (sum(isnan(epcs_all(i,:))) == 0);
            end
            
            epcs_all_cell(:) = (mat2cell(epcs_all, ones(size(epcs_all,1),1), size(epcs_all,2)))';

        
        miniepochs{ch}.data =  epcs_all_cell;
        miniepochs{ch}.miniepoch_labels = distance2sleep_onset;
        miniepochs{ch}.notnan = miniepoch_notnan;
        miniepochs{ch}.idxnotnan = find(miniepoch_notnan == 1); 
        miniepochs{ch}.sumnotnan = sum(miniepoch_notnan,2);
        miniepochs{ch}.isart = isart;
        miniepochs{ch}.miniepochs = epcs_all_cell;
        
        
        end
    
    end 
    epoched_table.EEGData = miniepochs;

    sbj2.sleep_onset.(col_name){k} = epoched_table;


end 

end


