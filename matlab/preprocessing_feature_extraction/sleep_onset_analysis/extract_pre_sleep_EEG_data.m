function presleep_EEG_mat = extract_pre_sleep_EEG_data(sbj2, duration)
% Extract the pre-sleep-onset EEG signal from the Patient class, and output
% it as a matrix (num_chs x num_samps), where num_chs -> number of EEG
% channels, num_samps -> number of datapoints in EEG period
% (period_duration_in_minutes x 60 seconds x sampling rate).
%
%
% Author: Anastasia Ilina
%
%% Function inputs:
%
% sbj2:                         the patient object for processing
% table_idx:                    the indexes from EEG_table to extract the
% labels:                       lables for the samples 
% 
%
%% Log of code:
%
% 29/06/2023 - Created by Anastasia Ilina 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if subject discarded
if sbj2.Sbj_discard
    disp('Subject discarded, no samples returned')
    presleep_EEG_mat = [];
    return
end

if isempty(sbj2.sleep_onset)
    error('Pre-sleep-onset data is not present in the subject, it might not have been extracted yet. Stop Analysis')
end 


sbj_id = sbj2.patient_original_id; 
Fs = sbj2.EEG_sampling_freq;
max_epochs = Fs * 60 * duration;



%sleep_scoring = sbj2.sleep_onset.PreSleepOnsetPeriodData{1,1}.SleepScoring{1,1};


%variable_names = {'Sbj_ID', 'Age', 'ifCleanOnset','SleepStaging', 'Channel_Name', 'EEG_data'};

%var_types = [];
%for var_type = 1:length(variable_names)
%    if var_type == length(variable_names) || var_type == length(variable_names) - 2
 %       var_types = [var_types, "cell"];
 %   elseif var_type == length(variable_names) - 1
  %      var_types = [var_types, "string"];
  %  else
  %      var_types = [var_types, "double"];
  %  end
%end 

%channel_table = table('Size', [num_chs, length(variable_names)], 'VariableTypes', var_types, 'VariableNames',variable_names);

%samp_table.Properties.VariableNames = variable_names;

%for samp_idx  = 1:num_samples
%    samp_now = samp_sbj_final_r{samp_idx};
%   samp_table.Label(samp_idx) = samp_now.Main_label;
%   samp_table.Sbj_ID(samp_idx) = samp_now.Sbj_id;
%    samp_table.Age(samp_idx) = sbj2.age;
%   samp_table.SleepStage(samp_idx) = sleep_scoring_samples(samp_idx);
%    samp_table.ifCleanOnset(samp_idx) = sbj2.sleep_onset.ifCleanOnset(1);
%    samp_table(samp_idx, 6:end) = array2table(samp_now.ft_vec); 
% end


presleep_EEG_mat = [];

EEGData = sbj2.sleep_onset.PreSleepOnsetPeriodData{1,1}.EEGData{1,1};
num_chs = size(EEGData, 1);

for ch = 1:num_chs
    if length(EEGData.EEG_data{ch}) == max_epochs
        presleep_EEG_mat(1, ch, :) = EEGData.EEG_data{ch};
    else
        %disp(size(EEGData.EEG_data{ch}))
        extra_padding = NaN([max_epochs - length(EEGData.EEG_data{ch}), 1]);
        %disp(size(extra_padding))
        padded_eeg = [extra_padding; EEGData.EEG_data{ch}];
        if length(padded_eeg) ~= max_epochs
            error('Padded EEG is not the right size')
        end 
        presleep_EEG_mat(1, ch, :) = padded_eeg; 
    end 
end 

end 



