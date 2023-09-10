function sbj2 = prepare_EEG_for_higher_order_fts(sbj2, table_idx, labels)
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
%  sbj2:                         the patient object for processing
% 
%
%% Log of code:
%
% 11/06/2023 - Created by Anastasia Ilina
% 16/06/2023 - Finished by Anastasia Iline
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check if subject discarded
if sbj2.Sbj_discard
    disp('Subject discarded, no samples returned')
    
    return
end

num_labels = length(labels);
if length(table_idx) ~= num_labels
    error('The number of labels should correspond to each table index')
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
    sbj2.Sbj_discard = 1;
    return 
end

% Check if the channels have the same number of samples

for lb = 1:num_labels
    tb_idx_now = table_idx{lb};
    num_samp_1 = length(sbj2.EEG_table.(tb_idx_now){1});
    for ch = 1: num_chs
        num_samp_ch = length(sbj2.EEG_table.(tb_idx_now){ch});
        if num_samp_ch ~= num_samp_1
            disp('Number of samples between channels is different. Will discard the participant')
            sbj2.Sbj_discard = 1; 
            return 
        end
    end
end 
           
            
% Create a table to work with (containing all the samples (num_ch x
% sample_lenght) and their class (Awake/Sleep_Onset)
vartypes = [];
for x = 1: num_samp_1
    vartypes = [vartypes; "double"]; 
end 
%EEG_sample_table = table('Size',[num_labels, num_samp_1], 'VariableTypes', vartypes , 'RowNames', labels);
EEG_sample_cell = cell([num_labels, num_samp_1]);


for lb = 1:num_labels
    
    tb_idx_now = table_idx{lb};  %iterating through the table indices and labels input to the code e.g. Awake -> Sleep Onset -> aslee

    label_now = labels{lb};    
    
    
        for i = 1:num_samp_1 %iterate for all of the samples for that type of EEG
            
            EEG_sample_all_chs = [];
            
            for ch = 1:num_chs
        
             % Check channel discard
             %if sbj2.EEG_table.discarded(ch)
                %    continue
                %end

                current_sample_EEG = sbj2.EEG_table.(tb_idx_now){ch}{i}; %take out the number of EEG traces (samples of 4s long) corresponding to the table index for that patient index
               
                
                % see how many samples there are for that type of EEG traces for that patient
                  % Clear object
                %num_sample = num_sample + 1;
                EEG_sample_all_chs(ch, :) = current_sample_EEG;
                
            
            end 
            sample_len = max(size(current_sample_EEG)); 
           
            %EEG_sample_all_chs_table = array2table(EEG_sample_all_chs); 
            %EEG_sample_all_chs_table.Properties.VariableNames = 1:sample_len;
            %EEG_sample_all_chs_table.Properties.RowNames = channel_names;
            EEG_sample_cell{lb,i} = EEG_sample_all_chs;
            EEG_sample_table = cell2table(EEG_sample_cell);
            %EEG_sample_table.Properties.VariableNames = 1:num_samp_1;
            EEG_sample_table.Properties.RowNames = labels;
            
          
        end 
end 

sbj2.higher_order_EEG_table = EEG_sample_table;