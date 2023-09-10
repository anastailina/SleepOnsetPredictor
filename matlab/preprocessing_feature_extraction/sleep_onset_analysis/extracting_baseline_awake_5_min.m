%% Test of patient class
cd '/home/anastasia/Dropbox/missing data'

load MESA_ID.mat

% Load demographic information
load MESA_Demo.mat
load ECG_quality.mat
load EEG_quality.mat
load 'Sbjs_filtered_time (2).mat'

EEG_headers = {'EEG1','EEG2','EEG3'};  % The labels for the 3 EEG channels
EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};

ECG_header = {'EKG'};  % ECG channel label

%load('/home/anastasia/filter for participants/Sbjs_filtered_time (2).mat')
cleanMESA_ID = MESA_ID(sbj_remained);
cleanMESA_Age = MESA_Age(sbj_remained);
cleanMESA_Race = MESA_Race(sbj_remained);
cleanMESA_Gender = MESA_Gender(sbj_remained);
num_sbjs = length(cleanMESA_ID);    

cleanEEG_FzCz_qual = EEG_FzCz_qual(sbj_remained);
cleanEEG_CzOz_qual = EEG_CzOz_qual(sbj_remained);
cleanEEG_C4M1_qual = EEG_C4M1_qual(sbj_remained);
cleanECG_qual = ECG_qual(sbj_remained);

sample_all = {};
timelimit = 8;

for sbj_id = clean_sbj_IDs%:length(cleanMESA_ID)
    cd /mnt/LongTermStorage/mesa/polysomnography/edfs
 
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
    sbj2 = cardio_art_rm(sbj2,sbj2.ifcard_marker,'Artefact_removed_signal');
end 

labels = {'Baseline_Awake'};

for sbj_id = 1:length(cleanMESA_ID)

    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/sbj'
    
    sbj_file = sbj_filename(cleanMESA_ID(sbj_id));
    
    load(sbj_file)
 
    sbj2 = sleep_staging(sbj2,2,'N2', 'Baseline Awake Epochs', 10); % each epoch is 30s, want to extract 5 minutes

    sbj2 = baseline_awake_extract(sbj2,0);

    sbj2 = sleep_analysis(sbj2,30,36000);

    sbj2 = ft_epoch_ext(sbj2,6,'Baseline_Awake_Period',30, 'Overlap', 0);
    
    if sbj2.Sbj_discard
        continue
    end 

    sbj2 = prepare_EEG_for_higher_order_fts(sbj2, {'Baseline_Awake_Period_epochs'}, labels);
    
    if sbj2.Sbj_discard
        continue
    end 
    %sbj2 = ft_epoch_ext(sbj2,6,'Onset_period',4,'Overlap',0.5);
    
    %% Save the resulting patient files 

    % change the path to a desired folder for saving 
    %cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    %filename = sbj_filename(cleanMESA_ID(sbj_id));
    %save (filename, 'sbj2')

%end


 %% Feature extraction 
%for sbj_id = clean_sbj_IDs % 1:length(cleanMESA_ID)
    %cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    %filename = sbj_filename(cleanMESA_ID(sbj_id));
    %load (filename)

    % Notice that samp_sbj2 will be a cell structure containing all the samples
    % from that subject

    samp_sbj2 = sample_gen(sbj2,{'Baseline_Awake_Period_epochs'},{'Baseline_Awake'},0);
    

    % In real analysis where multiple subjects are available, you can combine
    % them into one sample cell strucutre in order to compute features

    %! Remember to compile the Catch-22 files!
    %run('mexAll.m')
    
    % Extract features per channel 
    
    samp_sbj2 = ft_eval_sleep_onset(samp_sbj2,'new_sleep_EEG_ft_dscrp.txt',0);
    %samp_sbj2 = ft_eval_sleep_onset_without_gamma(samp_sbj2,'new_sleep_EEG_ft_dscrp_without_gamma.txt',0);

    % Store all samples in one variable (uncomment if needed)
    
    %  sample_all = [samp_sbj2, samp_sbj2];
    

    %% Save the extracted features
    
    % change the path to a desired folder for saving 
    
    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/unichannel'
    
    filename = sample_filename(cleanMESA_ID(sbj_id));
    
    save (filename, 'samp_sbj2')

%end 


%% Higher-order Feature extraction 
%for sbj_id =clean_sbj_IDs %1:length(cleanMESA_ID)
    
    %cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    %filename = sbj_filename(cleanMESA_ID(sbj_id));
    %load (filename)

    % Extract multichannel features

    

    high_order_fts_sbj2 = extract_high_order_features(sbj2, {'Baseline_Awake_Period_epochs'},labels);

    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/multichannel'
    
    ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));
    
    save (ho_filename, 'high_order_fts_sbj2')
    

end


%% Join together unichannel and multichannel features 
discarded_sbjs = [];

% Load the files
for sbj_id =1:length(cleanMESA_ID)


    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/unichannel' 
    
    filename = sample_filename(cleanMESA_ID(sbj_id));

    if isfile(filename)
        load (filename)
    else
        disp(['No unichannel sample file for Sbj. ', num2str(cleanMESA_ID(sbj_id))])
        discarded_sbjs = [discarded_sbjs,sbj_id];
        continue
    end 

    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/multichannel'  
    
    ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));

    if isfile(ho_filename)
        load (ho_filename)
    else
        disp(['No unichannel sample file for Sbj. ', num2str(cleanMESA_ID(sbj_id))])
        discarded_sbjs = [discarded_sbjs,sbj_id];
        continue
    end 

    if isempty(high_order_fts_sbj2) || isempty(samp_sbj2)
        discarded_sbjs = [discarded_sbjs,sbj_id];
        disp(['Sbj ', num2str(cleanMESA_ID(sbj_id)), ' is discarded as one of the sample cells is empty'])

        continue
    end 
    
    EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};


    % If we are choosing to average sample features across 3 channels, set
    % ifaveraged to 1

    ifaveraged = 1;

    if ifaveraged
        samp_sbj2_updated = average_features_across_chs(samp_sbj2, 3);
       
    else
        samp_sbj2_updated = create_unique_features(samp_sbj2, EEG_chs);
    end 
    
    % Join the samples together

    samp_sbj2_final = join_fts_together(samp_sbj2_updated, high_order_fts_sbj2, ifaveraged);

    
    % Save the final samples
    
    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/final'
    
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    
    save (final_filename, 'samp_sbj2_final')

    disp(['Samples from Sbj ', num2str(cleanMESA_ID(sbj_id)), ' have been joined'])

end 


%% Join them into a table


for sbj_id = 1:length(cleanMESA_ID)

    %cd  /mnt/FastData/anastasia/sleep_onset/samp_sbj_1min_proof_of_concept/final_samples

    %cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/samp_sbj2/final_each_channel'
    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/final'
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    

    if isfile(final_filename)
        load (final_filename)
    else
        disp(['No sample file for Sbj. ', num2str(cleanMESA_ID(sbj_id))])
        continue
    end 

    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/sbj'

    sbj2_filename = sbj_filename(cleanMESA_ID(sbj_id));

    if isfile(sbj2_filename)
        load (sbj2_filename)
    else 
        continue 
    end 

    
    % create a sample table

    sbj_sample_table = create_sample_table(samp_sbj2_final, sbj2, 6); 


    % join tables across the participants into one big table 

    if  sbj_id == 2 
        all_sbj_sample_table = sbj_sample_table;
    else
        all_sbj_sample_table =[all_sbj_sample_table; sbj_sample_table];
    end 
    disp(['Joining completed for Sbj. ', num2str(cleanMESA_ID(sbj_id)), '. Saving the file.'])
end 


 % export the table to csv format (save as a csv file);
cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s'
%writetable(all_sbj_sample_table, 'regression_so_samples_unique_ch_samples.csv')
writetable(all_sbj_sample_table, 'baseline_awake_5min_6sec.csv')

cd /home/anastasia/Dropbox
writetable(all_sbj_sample_table, 'baseline_awake_5min_6sec.csv')



%% Join all of the samples together 

%sbj_IDs = 1:length(cleanMESA_ID);
%clean_sbj_IDs = sbj_IDs(~ismember(sbj_IDs, discarded_sbjs));

sample_all = {};

for sbj_id = clean_sbj_IDs

    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/baseline_awake_5min_6s/final'
    

    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));

    if isfile(final_filename)
        load(final_filename)
    else
        disp(['No final sample for Sbj ', num2str(cleanMESA_ID(sbj_id))])
        continue
    end 

    sample_all = [sample_all, samp_sbj2_final];

    disp(['Samples from Sbj ', num2str(cleanMESA_ID(sbj_id)), ' concatenated'])

end 


