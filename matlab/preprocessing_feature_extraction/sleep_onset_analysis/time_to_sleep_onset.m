%% Test of patient class

%cd '/Users/anastasia/Dropbox'/'missing data'/
cd 'home/anastasia/Dropbox/missing data'

load MESA_ID.mat

% Load demographic information
load MESA_Demo.mat
load ECG_quality.mat
load EEG_quality.mat
%load 'Sbjs_filtered_time (2).mat'
%load sbj_final.mat
load('Users/anastasia/Dropbox/My code/new pipeline/filter for participants/sbj_final.mat', 'sbj_final')
load '/home/anastasia/Dropbox/filter for participants/Sbjs_filtered_time (4).mat'


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

for sbj_id = 1:n %length(cleanMESA_ID)

    cd '/Users/anastasia/Dropbox/sbj examples/raw_sbjs'
    
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    load(filename)

    sbj2 = cardio_art_rm(sbj2,sbj2.ifcard_marker,'Artefact_removed_signal');

    sbj2 = sleep_staging(sbj2,2,'N2');

    sbj2 = period_extract(sbj2,1);

    sbj2 = sleep_analysis(sbj2,10,36000);
    
    sbj2 = identify_sleep_onset(sbj2,2,'N2'); % identify the start of the sleep period and bedtime 

    sbj2 = extract_presleep_onset_EEG(sbj2, 30); % extract continious pre-sleep-onset EEG period 
    
    sbj2 = pre_sleep_onset_EEG_epoching(sbj2, 6, 0, 30); % break it up into 6 second epochs 


    %% Save the resulting patient files 

    % change the path to a desired folder for saving 
    cd '/Users/anastasia/Dropbox/sbj examples/sbj2'
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    save (filename, 'sbj2')

end

%% Feature extraction (Single channel features)

for sbj_id = 1:n %length(cleanMESA_ID)

    %cd /mnt/FastData/anastasia/sleep_onset/time_to_sleep_regression/30_mins/sbj2
    cd '/Users/anastasia/Dropbox/sbj examples/sbj2'
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    load (filename)
    
    EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};
    ifcardiac = 0;
    samp_sbj2 = sample_gen_sleepons_regression(sbj2, EEG_chs, ifcardiac);
    samp_sbj2 = ft_eval_sleep_onset(samp_sbj2,'new_sleep_EEG_ft_dscrp.txt',0);

    %% Save the resulting sample files 
     % change the path to a desired folder for saving 
    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/unichannel'
    %cd
    %/mnt/FastData/anastasia/sleep_onset/time_to_sleep_regression/30_mins/samp_sbj2/
    filename = sample_filename(cleanMESA_ID(sbj_id));
    save (filename, 'samp_sbj2')
    
end 

%% Feature extraction (multichannel features)


for sbj_id = 1: n %length(cleanMESA_ID)
    
    %cd /mnt/FastData/anastasia/sleep_onset/time_to_sleep_regression/30_mins/samp_sbj2
    cd '/Users/anastasia/Dropbox/sbj examples/sbj2'

    filename = sbj_filename(cleanMESA_ID(sbj_id));

    load (filename)

    sbj2 = prepare_EEG_for_higher_order_fts_regression(sbj2, EEG_chs);


    high_order_fts_sbj2 = extract_ho_fts_so_regression(sbj2,1);
    
    if ~isempty(high_order_fts_sbj2)
        cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/high_order'

        ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));
        save (ho_filename, 'high_order_fts_sbj2')
    end


end 

%% Join together unichannel and multichannel features 

% Load the files
for sbj_id = 1:n%length(cleanMESA_ID)

    %cd  /mnt/FastData/anastasia/sleep_onset/samp_sbj_1min_proof_of_concept/unichannel_ft 
    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/unichannel'

    filename = sample_filename(cleanMESA_ID(sbj_id));
    if isfile(filename)
        load (filename)
    else
        continue
    end 


   % cd  /mnt/FastData/anastasia/sleep_onset/samp_sbj_1min_proof_of_concept/multichannel_ft 
    
    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/high_order'

    ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));


    if isfile(ho_filename)
        load (ho_filename)
    else
        continue
    end 

    
    EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};


    % If we are choosing to average sample features across 3 channels, set
    % ifaveraged to 1

    ifaveraged = 1;

    if ifaveraged
        samp_sbj2_updated_r = average_features_across_chs_regression(samp_sbj2, 3);
       
    else
        samp_sbj2_updated_r = create_unique_features(samp_sbj2, EEG_chs);
    end 
    
    
    % Join the samples together

    samp_sbj2_final_r = join_fts_together(samp_sbj2_updated_r, high_order_fts_sbj2);

    
    % Save the final samples
    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/final_samples'
    
   % cd  /mnt/FastData/anastasia/sleep_onset/samp_sbj_1min_proof_of_concept/final_samples
    
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    
    save (final_filename, 'samp_sbj2_final_r')

end 

%% Prepare the samples in format suitable for python

for sbj_id = 1:n+1 %length(cleanMESA_ID)

    %cd  /mnt/FastData/anastasia/sleep_onset/samp_sbj_1min_proof_of_concept/final_samples

    cd '/Users/anastasia/Dropbox/sbj examples/samp_sbj2/final_samples'
    
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    

    if isfile(final_filename)
        load (final_filename)
    else
        continue
    end 

    cd '/Users/anastasia/Dropbox/sbj examples/sbj2'

    sbj2_filename = sbj_filename(cleanMESA_ID(sbj_id));

    if isfile(sbj2_filename)
        load (sbj2_filename)
    else 
        continue 
    end 

    
    % create a sample table

    sbj_sample_table = create_sample_table(samp_sbj2_final_r, sbj2, 6); 


    % join tables across the participants into one big table 

    if  sbj_id == 1 
        all_sbj_sample_table = sbj_sample_table;
    else
        all_sbj_sample_table =[all_sbj_sample_table; sbj_sample_table];
    end 

end 


 % export the table to csv format (save as a csv file);
 cd /Users/anastasia/Dropbox/Sleeper/sleep_onset_analysis
 writetable(all_sbj_sample_table, 'regression_so_samples_61_subjects.csv')

 

