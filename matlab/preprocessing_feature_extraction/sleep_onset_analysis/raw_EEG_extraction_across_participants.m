%% Test of patient class

cd '/Users/anastasia/Dropbox/missing data'
%cd 'home/anastasia/Dropbox/missing data'

load MESA_ID.mat

% Load demographic information
load MESA_Demo.mat
load ECG_quality.mat
load EEG_quality.mat
%load 'Sbjs_filtered_time (2).mat'
%load sbj_final.mat
load('Users/anastasia/Dropbox/My code/new pipeline/filter for participants/sbj_final.mat', 'sbj_final')
load 'Users/anastasia/Dropbox/filter for participants/Sbjs_filtered_time (4).mat'


EEG_headers = {'EEG1','EEG2','EEG3'};    % The lables for the 3 EEG channels
EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};

ECG_header = {'EKG'};               % ECG channel label

EEG_quality_sbj = zeros(length(EEG_chs),1);
EEG_quality_sbj(1) = EEG_FzCz_qual(sbj_id);
EEG_quality_sbj(2) = EEG_CzOz_qual(sbj_id);
EEG_quality_sbj(3) = EEG_C4M1_qual(sbj_id);

ECG_quality_sbj = ECG_qual(sbj_id);

cleanMESA_ID = MESA_ID(sbj_final);
cleanMESA_Age = MESA_Age(sbj_final);
cleanMESA_Race = MESA_Race(sbj_final);
cleanMESA_Gender = MESA_Gender(sbj_final);
num_sbjs = length(cleanMESA_ID);    

cleanEEG_FzCz_qual = EEG_FzCz_qual(sbj_final);
cleanEEG_CzOz_qual = EEG_CzOz_qual(sbj_final);
cleanEEG_C4M1_qual = EEG_C4M1_qual(sbj_final);
cleanECG_qual = ECG_qual(sbj_final);

timelimit = 8;
sample_all = {};

for sbj_id = 1:length(cleanMESA_ID)
 
    EEG_quality_sbj = zeros(length(EEG_chs),1);
    EEG_quality_sbj(1) = cleanEEG_FzCz_qual(sbj_id);
    EEG_quality_sbj(2) = cleanEEG_CzOz_qual(sbj_id);
    EEG_quality_sbj(3) = cleanEEG_C4M1_qual(sbj_id);

    ECG_quality_sbj = cleanECG_qual(sbj_id);

    sbj1 = Patient('patient_original_id', cleanMESA_ID(sbj_id),...
                'age',cleanMESA_Age(sbj_id),...
                'race',cleanMESA_Race(sbj_id),...
                'gender',cleanMESA_Gender(sbj_id),...
                'file_prefix','mesa-sleep');

    sbj1 = sbj1.init();


    sbj1 = loadPSG(sbj1,'EEG_channel_names',EEG_chs,...
               'EEG_header_index',EEG_headers,...
               'ECG_header_index',ECG_header,...
               'Time_to_take',timelimit, 'EEG_quality',EEG_quality_sbj,...
               'ECG_quality',ECG_quality_sbj);

    sbj1 = readAnnot(sbj1); %do not fully understand what this function does 

    sbj1 = readRmarker(sbj1,'Epoch_idx','epoch',...
                   'Rpoint_idx','RPoint',...
                   'QRS_start_idx','Start',...
                   'QRS_end_idx','End');

      %% Pre-processing

    sbj1 = quality_check(sbj1,3);

    sbj1 = general_art_rm(sbj1,'NormalRMS',2,'Signal');

    sbj2 = sbj1;
end 

n = 61;
all_data = [];

for sbj_id = 1:n %length(cleanMESA_ID)
    
    presleep_eeg_data = [];

    cd '/Users/anastasia/Dropbox/sbj examples/raw_sbjs'
    
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    load(filename)

    sbj2 = cardio_art_rm(sbj2,sbj2.ifcard_marker,'Artefact_removed_signal');

    sbj2 = sleep_staging(sbj2,2,'N2');

    sbj2 = period_extract(sbj2,1);

    sbj2 = sleep_analysis(sbj2,10,36000);
    
    sbj2 = identify_sleep_onset(sbj2,2,'N2'); % identify the start of the sleep period and bedtime 

    sbj2 = extract_presleep_onset_EEG(sbj2, 30); % extract continious pre-sleep-onset EEG period 
    
    presleep_eeg_data = extract_pre_sleep_EEG_data(sbj2, 30);
   
    if ~isempty(presleep_eeg_data)
        
        all_data = cat(1, all_data, presleep_eeg_data);

        disp(['EEG data from Sbj ', num2str(sbj_id), ' out of ', num2str(n),' was concatenated'])
    end
end

cd '/Users/anastasia/Dropbox/sbj examples'
save ('allEEGdata.mat','all_data')
