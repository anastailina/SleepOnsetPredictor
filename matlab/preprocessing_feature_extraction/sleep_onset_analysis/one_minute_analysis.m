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
%end 

%for sbj_id = 1:length(cleanMESA_ID)

    %cd /mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/sbj
    
    %filename = sbj_filename(cleanMESA_ID(sbj_id));
 
    sbj2 = sleep_staging(sbj2,2,'N2');

    sbj2 = period_extract(sbj2,1);

    sbj2 = sleep_analysis(sbj2,10,36000);

    sbj2 = ft_epoch_ext(sbj2,5,'Awake_period',4,'Overlap',0.5);
    sbj2 = ft_epoch_ext(sbj2,5,'Onset_period',4,'Overlap',0.5);
    
    %% Save the resulting patient files 

    % change the path to a desired folder for saving 
    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    save (filename, 'sbj2')

end


 %% Feature extraction 
for sbj_id = clean_sbj_IDs % 1:length(cleanMESA_ID)
    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    load (filename)

    % Notice that samp_sbj2 will be a cell structure containing all the samples
    % from that subject

    samp_sbj2 = sample_gen(sbj2,{'Awake_period_epochs','Onset_period_epochs'},{'Awake','Onset'},1);
    

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
    
    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/unichannel'
    
    filename = sample_filename(cleanMESA_ID(sbj_id));
    
    save (filename, 'samp_sbj2')

end 


%% Higher-order Feature extraction 
for sbj_id =clean_sbj_IDs %1:length(cleanMESA_ID)
    
    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/sbj'
    filename = sbj_filename(cleanMESA_ID(sbj_id));
    load (filename)

    % Extract multichannel features

    labels = {'Awake','Onset'};

    sbj2 = prepare_EEG_for_higher_order_fts(sbj2, {'Awake_period_epochs','Onset_period_epochs'},labels);

    high_order_fts_sbj2 = extract_high_order_features(sbj2, {'Awake_period_epochs','Onset_period_epochs'},labels);

    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/multichannel'
    
    ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));
    
    save (ho_filename, 'high_order_fts_sbj2')
    

end


%% Join together unichannel and multichannel features 
discarded_sbjs = [];

% Load the files
for sbj_id = clean_sbj_IDs %1:length(cleanMESA_ID)

    cd '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/unichannel' 
    
    filename = sample_filename(cleanMESA_ID(sbj_id));
    
    load (filename)


    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/multichannel'  
    
    ho_filename = higher_order_sample_filename(cleanMESA_ID(sbj_id));
    
    load (ho_filename)
    
    
   
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

    samp_sbj2_final = join_fts_together(samp_sbj2_updated, high_order_fts_sbj2);

    
    % Save the final samples
    
    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/final'
    
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    
    save (final_filename, 'samp_sbj2_final')

    disp(['Samples from Sbj ', num2str(cleanMESA_ID(sbj_id)), ' have been joined'])

end 

%% Join all of the samples together 

%sbj_IDs = 1:length(cleanMESA_ID);
%clean_sbj_IDs = sbj_IDs(~ismember(sbj_IDs, discarded_sbjs));

sample_all = {};

for sbj_id = clean_sbj_IDs

    cd  '/mnt/LongTermStorage/anastasia/mesa_analysed/thesis 2023/sleep onset/proof_of_concept/samp_sbj2/final'
    
    final_filename = final_sample_filename(cleanMESA_ID(sbj_id));
    
    load(final_filename)

    sample_all = [sample_all, samp_sbj2_final];

    disp(['Samples from Sbj ', num2str(cleanMESA_ID(sbj_id)), ' concatenated'])

end 





%Spliting the sample_all into 8 parts

%split = 8;
num_smpl = size(sample_all,2);
%idx_split = floor(num_smpl/split);
%sample1 = sample_all(1:idx_split);
%%sample3 = sample_all(2*idx_split+1: idx_split*3);
%sample4 = sample_all(3*idx_split+1: idx_split*4);
%sample5 = sample_all(4*idx_split+1: idx_split*5);
%sample6 = sample_all(5*idx_split+1: idx_split*6);
%sample7 = sample_all(6*idx_split+1: idx_split*7);
%sample8 = sample_all(7*idx_split+1: end);


%Extracting only features form the sample_all

feature_sbj_all = {}
for i = 1: num_smpl
    feature_subj_all{i}.ft_vec = sample_all{i}.ft_vec;
    feature_subj_all{i}.ft_dscrp = sample_all{i}.ft_dscrp;
    feature_subj_all{i}.Main_label = sample_all{i}.Main_label;
    feature_subj_all{i}.Channel_label = sample_all{i}.Channel_label;
    feature_subj_all{i}.Sbj_id = sample_all{i}.Sbj_id;
    feature_subj_all{i}.Fs = sample_all{i}.Fs;
end



%% Basic machine learning

% create a matrix of all features and the corresponding label array
num_obs = size(feature_subj_all,2);
num_ft = size(feature_subj_all{1}.ft_vec,2);

labels = {};

feature_space = [];
for smp = 1:num_obs
    feature_space(smp,:) = feature_subj_all{smp}.ft_vec; %creating a variable containing the feature space
    labels{smp} = feature_subj_all{smp}.Main_label;
end

cat_labels = categorical((labels)');

save ('feature_space_without_gamma.mat', 'feature_space')
save ('labels_without_gamma.mat', "labels")
save ('feature_subj_all_without_gamma.mat','feature_subj_all')
 
labels_test = [];
labels_train = [];
feature_space_train = [];
feature_space_test = [];
% Splitiing it up into test and training data
p = 0.8; %percentage of the data for training
rand_idx = randperm(num_obs);
idx_train = rand_idx(1:round(p*num_obs));
idx_test = rand_idx((round(p*num_obs)+1):end);
feature_space_train = feature_space(rand_idx(idx_train),:);
labels_train = labels(:,rand_idx(idx_train));

feature_space_test = feature_space(rand_idx(idx_test),:);
labels_test = labels(:,rand_idx(idx_test));

%fit the SVM classifier (train)
SVMModel = fitcsvm(feature_space_train,labels_train);
CVSVMModel = crossval(SVMModel); %Cross-validate the SVM classifier. By default uses 10-fold cross-validation
%Classifying new data with an SVM classifier
[pred_label,score] = predict(SVMModel,feature_space_test);


for lbl = 1:size(labels_test,2)
    pred_label_str(lbl) = string(pred_label{lbl});
end 

%Calculating the accuracy for linear SVM 
test_accuracy_for_iter = (sum(pred_label_str == labels_test)/length(labels_test))*100

%Extracting the top-performing features and performing hierarchical
%clusterig of the features
ft_params = feature_parameters(feature_subj_all);
cfnParams = GiveMeDefaultParams(labels);
[ifeat,testStat] = topFeatureHierClust(feature_space,cfnParams,ft_params, 40);
