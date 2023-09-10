function samp_table =create_sample_table(samp_sbj_final_r, sbj2, epoch_dur)
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


num_samples = length(samp_sbj_final_r);
sbj_id = samp_sbj_final_r{1,1}.Sbj_id;
scoring_epoch_duration = sbj2.Scoring_Epoch;
sleep_scoring = sbj2.sleep_onset.PreSleepOnsetPeriodData{1,1}.SleepScoring{1,1};
sleep_scoring_samples = [];

for i = 1: num_samples
    new_idx = i/(scoring_epoch_duration/epoch_dur);
    sleep_scoring_samples(i) = sleep_scoring(ceil(new_idx)); 
end 

variable_names = {'Label', 'Sbj_ID', 'Age', 'ifCleanOnset','SleepStage'};

variable_names = [variable_names, samp_sbj_final_r{1,1}.ft_dscrp];

var_types = [];
for var_type = 1:length(variable_names)
    var_types = [var_types, "double"];
end 

samp_table = table('Size', [num_samples, length(variable_names)], 'VariableTypes', var_types, 'VariableNames',variable_names);

samp_table.Properties.VariableNames = variable_names;

for samp_idx  = 1:num_samples
    samp_now = samp_sbj_final_r{samp_idx};
    samp_table.Label(samp_idx) = 1;
    samp_table.Sbj_ID(samp_idx) = samp_now.Sbj_id;
    samp_table.Age(samp_idx) = sbj2.age;
    samp_table.SleepStage(samp_idx) = sleep_scoring_samples(samp_idx);
    samp_table.ifCleanOnset(samp_idx) = sbj2.sleep_onset.ifCleanOnset(1);
    samp_table(samp_idx, 6:end) = array2table(samp_now.ft_vec); 
end


