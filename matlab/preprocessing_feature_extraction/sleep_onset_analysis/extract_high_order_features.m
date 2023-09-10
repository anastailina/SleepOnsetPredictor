function  sampleObj_cell = extract_high_order_features(sbj2, table_idx, labels)

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
% table_idx:                    the indexes from EEG_table to extract the
% labels:                       lables for the samples 
% 
%
%% Log of code:
%
% 11/06/2023 - Created by Anastasia Ilina 
% 15/06/2023 - Finished by Anastasia Ilina 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check if subject discarded
if sbj2.Sbj_discard
    disp('Subject discarded, no samples returned')
    sampleObj_cell = {};
    return
end

num_labels = length(labels);
if length(table_idx) ~= num_labels
    error('The number of labels should correspond to each table index')
end

if isempty(sbj2.higher_order_EEG_table)
    disp('The higher order table is empty. Not proceeding, returning empty cell')
    sampleObj_cell = {};
    return
    
end 

num_chs = size(sbj2.EEG_table,1); % number of channels   
sbj_id = sbj2.patient_original_id;
eeg_Fs = sbj2.EEG_sampling_freq;
num_ft = 2;
% Check if none of the channels are discarded


EEG_sample_table = sbj2.higher_order_EEG_table;
sampleObj_cell = {};

sample_idx = 0;
for lb = 1:length(labels)
    label = labels{lb}; 

    num_samples = max(size(EEG_sample_table(lb,:))); 
   

    for i = 1:num_samples 
        
        samp_obj = {};
        sample_idx = sample_idx+1;

        EEG_sample = EEG_sample_table(lb, i); 

        EEG_matrix = table2array(EEG_sample);

        %[Oinfo,Sinfo,Red,Syn] = high_order(EEG_matrix{1,1},num_chs);

        % [oinfo,sinfo] = soinfo_from_covmat(covmat,T) -> instead of
        % cov(EEG_matrix{1,1}') use correlation matrix -> corrcoef
        [~,est_covmat] = data2gaussian(EEG_matrix{1,1}');
        [oinfo,sinfo] = soinfo_from_covmat(est_covmat,1280);
        
        samp_obj = Sample(sbj_id);
        samp_obj.sampID = i;
        samp_obj.data_EEG = EEG_matrix; %input teh 
        samp_obj.Fs = eeg_Fs; %input the sampling frequency in there
        samp_obj.Main_label = label; % input the label in there
        samp_obj.Channel_label = EEG_sample.Properties.RowNames;
        samp_obj.ft_vec = [oinfo,sinfo];
        samp_obj.ft_dscrp = {"O-information", "S-information"};
        
        ft = samp_obj.ft_vec;

        %Checking if the subject has the right snumber of the extracted
        %features 
        if length(ft) ~= num_ft
            samp_obj.ft_quality = 0;
            warning(['Sample No.',num2str(i),' from Subject No.',num2str(samp_obj.Sbj_id),' does not have the right feature number'])
            continue
        else
            %check if any extracted features of this subject contain NaN values
            samp_obj.ft_vec = ft;
            if sum(isnan(ft))>0
                samp_obj.ft_quality = 0;
                disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_obj.Sbj_id),' have NaN feature value'])
                continue
            end
            samp_obj.ft_quality = 1;

            %If all is okay,siognifiy it by the ft_quality variable = 1 and
            %print out that ft-eval went well 
            disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_obj.Sbj_id),' has been successfully calculated'])
        end
    
        sampleObj_cell{sample_idx} = samp_obj; %input the extracted features and corresponding iformation back into teh samples form the subject structure
    
    end    


end 






