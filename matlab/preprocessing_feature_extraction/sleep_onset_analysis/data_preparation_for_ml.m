cd '/Users/anastasia/Dropbox'/'missing data'/
%cd '/home/anastasia/Dropbox/missing data'

load MESA_ID.mat

% Load demographic information
load MESA_Demo.mat
load ECG_quality.mat
load EEG_quality.mat
load 'Sbjs_filtered_time (2).mat'
%load sbj_final.mat
%load('Users/anastasia/Dropbox/My code/new pipeline/filter for participants/sbj_final.mat', 'sbj_final')
%load '/home/anastasia/Dropbox/filter for participants/Sbjs_filtered_time (4).mat'


cd '/Users/anastasia/Dropbox/anastasia_analysis'
load sbj_remained_30mins_final.mat
load newly_discarded_sbjs.mat
sbj_ids = setdiff(sbjs_remained_30mins_final, newly_discarded_sbjs); 
load oriscore_asleep.mat
[~, indices] = ismember(sbj_remained, sbj_ids);
indices_30_min = find(indices);
oriscore_asleep_remained = oriscore_asleep(indices_30_min);



cleanMESA_ID = MESA_ID(sbj_ids);
cleanMESA_Age = MESA_Age(sbj_ids);
cleanMESA_Race = MESA_Race(sbj_ids);
cleanMESA_Gender = MESA_Gender(sbj_ids);
num_sbjs = length(sbj_ids);    


threshold_selected = 4;

scoring_epoch_duration = 30;
pre_sleep_duration_min = 30;
baseline_awake_duration_min = 5; 
epoch_duration_seconds = 6;
epoch_duration_minutes = epoch_duration_seconds / 60;

num_sbjs = length(sbj_ids);
for sbj_id = 1:num_sbjs
    
    load(['new_fts_subject_', num2str(sbj_ids(sbj_id)), '.mat']);
    load(['subject_', num2str(sbj_ids(sbj_id)), '.mat'])


    real_sbj_id = cleanMESA_ID(sbj_id);
    sbj_age = cleanMESA_Age(sbj_id);
    sbj_race= cleanMESA_Race(sbj_id);
    sbj_gender = cleanMESA_Gender(sbj_id);

    artifact_table = sbj_data_new.artifact_table;
    th_row = artifact_table(artifact_table.Threshold == threshold_selected, :);
    art_mask_pre_sleep_all = th_row.Artifact_Mask_Pre_Sleep_All{1,1};
    art_mask_awake_all = th_row.Artifact_Mask_Baseline_Awake_All{1,1};
    
    num_samples_presleep = length(art_mask_pre_sleep_all);
    num_samples_awake = length(art_mask_awake_all);

    sleep_scoring_all = oriscore_asleep_remained{sbj_id};
    [firstFiveScores, thirtyBeforeOnsetScores, ifcleanonset, time2sleep] = sleepOnsetAnalysis(sleep_scoring_all, pre_sleep_duration_min, baseline_awake_duration_min, scoring_epoch_duration);
    sleep_scoring_samples_pre_sleep = NaN(num_samples_presleep, 1);
    sleep_scoring_samples_awake = NaN(num_samples_awake, 1);

    % Broadcast sleep scores 
    for i = 1: num_samples_presleep
        new_idx = i/(scoring_epoch_duration/epoch_duration_seconds);
        sleep_scoring_samples_pre_sleep(i) = thirtyBeforeOnsetScores(ceil(new_idx)); 
    end

    for i = 1: num_samples_awake
        new_idx = i/(scoring_epoch_duration/epoch_duration_seconds);
        sleep_scoring_samples_awake(i) = firstFiveScores(ceil(new_idx)); 
    end 
    
    new_ft_mat_presleep = sbj_data_new.new_ft_mat_pre_sleep;
    old_ft_mat_presleep = sbj_data_new.old_fts_pre_sleep;
    ho_ft_mat_presleep = sbj_data_new.ho_ft_mat_pre_sleep;

    new_ft_mat_start = sbj_data_new.new_ft_mat_start;
    old_ft_mat_start = sbj_data_new.old_fts_start;
    ho_ft_mat_start = sbj_data_new.ho_ft_mat_start;

    old_ft_dscrp = [sbj_data.ft_dscpr(1:19), sbj_data.ft_dscpr(26:end)];
    new_ft_dscrp = sbj_data_new.ft_dscrp_new;
    all_ft_dscrp = [old_ft_dscrp, new_ft_dscrp];
    num_fts = length(all_ft_dscrp);

    ft_mat_sbj_presleep = NaN(num_fts, num_samples_presleep);
    ft_mat_sbj_awake = NaN(num_fts, num_samples_awake);

    for samp = 1:num_samples_presleep
        if art_mask_pre_sleep_all(samp) == 1
            continue 
        else 
            ft_vec_all_ch = [mean(old_ft_mat_presleep{samp}.ft_ch(:, 1:19), 'omitnan'), mean(old_ft_mat_presleep{samp}.ft_ch(:, 26:end), 'omitnan') mean(new_ft_mat_presleep{samp}.ft_ch, 'omitnan'), ho_ft_mat_presleep{samp}];
            ft_mat_sbj_presleep(:, samp) = ft_vec_all_ch;
        end 
    end 
    
    for samp = 1:num_samples_awake
        if art_mask_awake_all(samp) == 1
            continue 
        else 
            ft_vec_all_ch = [mean(old_ft_mat_start{samp}.ft_ch(:, 1:19), 'omitnan'), mean(old_ft_mat_start{samp}.ft_ch(:, 26:end), 'omitnan') mean(new_ft_mat_start{samp}.ft_ch, 'omitnan'), ho_ft_mat_start{samp}];
            ft_mat_sbj_awake(:, samp) = ft_vec_all_ch;
        end 
    end 

    
    % Create the labels array
    labels_pre_sleep = epoch_duration_minutes * (num_samples_presleep - 1:-1:0);
    labels_awake  = epoch_duration_minutes * (0:1:num_samples_awake-1);

    variable_names = {'Label', 'Sbj_ID', 'Age', 'Gender', 'Race', 'ifCleanOnset', 'Time2Sleep', 'SleepStage'};

    variable_names = [variable_names, all_ft_dscrp];

    var_types = [];
    for var_type = 1:length(variable_names)
        var_types = [var_types, "double"];
    end 



    sbj_table_presleep = table('Size', [num_samples_presleep, length(variable_names)], 'VariableTypes', var_types, 'VariableNames',variable_names);
    sbj_table_awake = table('Size', [num_samples_awake, length(variable_names)], 'VariableTypes', var_types, 'VariableNames',variable_names);

    sbj_table_presleep.Properties.VariableNames = variable_names;
    sbj_table_awake.Properties.VariableNames = variable_names;
    
    for samp_idx  = 1:num_samples_presleep
        sbj_table_presleep.Label(samp_idx) = labels_pre_sleep(samp_idx);
        sbj_table_presleep.Sbj_ID(samp_idx) = real_sbj_id;
        sbj_table_presleep.Age(samp_idx) = sbj_age;
        sbj_table_presleep.Gender(samp_idx) = sbj_gender;
        sbj_table_presleep.Race(samp_idx) = sbj_race;
        sbj_table_presleep.SleepStage(samp_idx) =  sleep_scoring_samples_pre_sleep(samp_idx);
        sbj_table_presleep.ifCleanOnset(samp_idx) = ifcleanonset;
        sbj_table_presleep.Time2Sleep(samp_idx) = time2sleep;
        sbj_table_presleep(samp_idx, 9:end) = array2table(ft_mat_sbj_presleep(:, samp_idx)'); 
    end

    for samp_idx  = 1:num_samples_awake
        sbj_table_awake.Label(samp_idx) = labels_awake(samp_idx);
        sbj_table_awake.Sbj_ID(samp_idx) = real_sbj_id;
        sbj_table_awake.Age(samp_idx) = sbj_age;
        sbj_table_awake.Gender(samp_idx) = sbj_gender;
        sbj_table_awake.Race(samp_idx) = sbj_race;
        sbj_table_awake.SleepStage(samp_idx) = sleep_scoring_samples_awake(samp_idx);
        sbj_table_awake.ifCleanOnset(samp_idx) = ifcleanonset;
        sbj_table_awake.Time2Sleep(samp_idx) = time2sleep;
        sbj_table_awake(samp_idx, 9:end) = array2table(ft_mat_sbj_awake(:, samp_idx)'); 
    end
    
    sbj_data_new.sbj_table_presleep = sbj_table_presleep;
    sbj_data_new.sbj_table_awake = sbj_table_awake;
    sbj_data_new.sleep_score_10_after_sleep_onset = sleep_scoring_all;
    sbj_data_new.time2sleep = time2sleep;
    sbj_data_new.real_sbj_id = real_sbj_id;
    sbj_data_new.sbj_gender = sbj_gender;
    sbj_data_new.sbj_race = sbj_race;

    save(['new_fts_subject_', num2str(sbj_ids(sbj_id)), '.mat'], 'sbj_data_new')


    % join tables across the participants into one big table 

    if  sbj_id == 2
        all_sbj_sample_table_pre_sleep = sbj_table_presleep;
        all_sbj_sample_table_awake = sbj_table_awake;
    else
        all_sbj_sample_table_pre_sleep =[all_sbj_sample_table_pre_sleep; sbj_table_presleep];
        all_sbj_sample_table_awake = [all_sbj_sample_table_awake; sbj_table_awake];
    end 
    disp(['Joining completed for Sbj. ', num2str(cleanMESA_ID(sbj_id))])
   
end  

%cd /home/anastasia/Dropbox
cd /Users/anastasia/Dropbox
writetable(all_sbj_sample_table_pre_sleep, 'renewed_data_presleep.csv')
writetable(all_sbj_sample_table_awake, 'renewed_data_awake.csv')


