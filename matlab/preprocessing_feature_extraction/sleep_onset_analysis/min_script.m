%% Preparations

% Subject group
clear all
cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All/anastasia_analysis'

load sbj_remained_30mins_final.mat
load ft_dscrp_new.mat
ft_dscrp_new_with_ho = [ft_dscrp_new, {'O-information', 'S-information'}];


sbjs_ids = sbjs_remained_30mins_final;

%% Continuous feature evaluation pipelines

% Original analysis parameters
Stage = 'N2';    % Stage onset to study
num_sbjs = length(sbjs_ids);

% % Set up Python
flag = int32(bitor(2, 8));
py.sys.setdlopenflags(flag);
%python_directory = '/home/anastasia/Python-3.10.11/python';
iflinux = 1;
%setup_python(python_directory, iflinux);
% 

% Continuous extraction parameters

epc_len = 6;            % Epoch length for feature evaluation
ovlap = 0;              % Overlap between continous epochs
post_sleep_per = 0;    % Post sleep period for feature evaluation
Fs = 256;
order = 2;
[b,a] = butter(order, [0.1,100]/(Fs/2), 'bandpass'); %pre-filtering via Butterworth bandpass
N = 100;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];
gamma_bd = [31, 100];
d_delta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',delta_bd(1),'CutoffFrequency2',delta_bd(2),'SampleRate',Fs);
d_theta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',theta_bd(1),'CutoffFrequency2',theta_bd(2),'SampleRate',Fs);
d_alpha = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',alpha_bd(1),'CutoffFrequency2',alpha_bd(2),'SampleRate',Fs);
d_beta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',beta_bd(1),'CutoffFrequency2',beta_bd(2),'SampleRate',Fs);
d_gamma = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',gamma_bd(1),'CutoffFrequency2',gamma_bd(2),'SampleRate',Fs);
iflinux = 1;
python_directory = '/home/anastasia/Python-3.10.11/python';

% Initialize ft_sbj_newfts_all and num_ftepcs_all_new outside the loop
ft_sbj_newfts_all = cell(num_sbjs, 1);
num_ftepcs_all_new = zeros(num_sbjs, 1);
stp_ftepcs_all_new = cell(num_sbjs, 1);
newly_discarded_sbjs = [];

for sbj_id = 377:477%num_sbjs
    %sbj_id = 1;


    disp(['Subject No.',num2str(sbjs_ids(sbj_id)),' started calculations'])
    % Start the timer
    tic;
    
    % load the subject
    load(['subject_', num2str(sbjs_ids(sbj_id)), '.mat'])
    

    % Feature evaluation
    % load the eeg data for the subject
    eeg_mat_now = sbj_data.eeg_mat;
    num_chs = max(size(eeg_mat_now));
    ifempty = zeros(1, num_chs);
    % Check whether all of the channels are present
    
    for ch_now = 1:num_chs
        ifempty(ch_now) = isempty(eeg_mat_now{ch_now}); 
    end 
    
    if any(ifempty)
        disp(['Subject No.',num2str(sbjs_ids(sbj_id)),'has less than 3 channels'])
        newly_discarded_sbjs = [newly_discarded_sbjs, sbjs_ids(sbj_id)];
        continue
    end  

    % load the starting points of the epoch from Junheng data
    pre_sleep_start_points = sbj_data.pre_sleep_start_points;
    pre_sleep_old_fts = sbj_data.pre_sleep_old_fts;
    start_recording_start_points = sbj_data.start_recording_start_points;
    start_recording_old_fts = sbj_data.start_recording_old_fts;


    % initialise counter for the channel index
    idx_realch = 0;

    % get the channel of our the eeg matrix    
       
    len_now = length(eeg_mat_now{1});
    epc_len_data = epc_len*Fs;
    ovlap_data = floor(epc_len_data*ovlap);
    rm_data = epc_len_data-ovlap_data;
        
    % calculate the total number of epochs given overlapping
    num_epcs_total_pre_sleep = length(pre_sleep_start_points);
    num_epcs_total_start_rec = length(start_recording_start_points);

    % initialise cell to store fts for extracted samples 
    ft_sbj_newfts_pre_sleep = cell(1, num_epcs_total_pre_sleep);
    ft_sbj_newfts_start_rec = cell(1, num_epcs_total_start_rec);
    

    % get eeg_indices array
    eeg_presleep_idx = pre_sleep_start_points(1):pre_sleep_start_points(num_epcs_total_pre_sleep) + epc_len_data - 1;
    eeg_start_rec_idx = start_recording_start_points(1):start_recording_start_points(num_epcs_total_start_rec)+epc_len_data - 1;
    
    % prealocate eeg_mats for higher order functions
    eeg_mat_pre_sleep = NaN(num_chs, length(eeg_presleep_idx)); 
    eeg_mat_start_rec = NaN(num_chs, length(eeg_start_rec_idx)); 

    % iterate through the channel
    for ch = 1:num_chs
        
        % increment the channel counter
        idx_realch = idx_realch+1;
        
        % get the channel of our the eeg matrix
        eeg_ck_pre_sleep = eeg_mat_now{idx_realch}(eeg_presleep_idx);
        eeg_ck_start = eeg_mat_now{idx_realch}(eeg_start_rec_idx);

        eeg_mat_pre_sleep(ch,:) = eeg_ck_pre_sleep; 
        eeg_mat_start_rec(ch,:) = eeg_ck_start;
        

         % increment the total number of the epochs that the subject has
        %num_ftepcs_all_new(sbj_id) = num_ftepcs_all_new(sbj_id) + num_epcs_total; 

        % Extract features for pre-sleep period
        [ft_sbj_newfts_ch_pre_sleep, stp_ftepcs_new_pre_sleep] = process_epochs(num_epcs_total_pre_sleep, ...
                                                                              rm_data, eeg_ck_pre_sleep, epc_len_data, ...
                                                                              pre_sleep_old_fts, ch, Fs, b, a, python_directory, ...
                                                                              d_delta, d_theta, d_alpha, d_beta, d_gamma, ...
                                                                              iflinux);
        
        
        % Extract features for  the first 5 minutes of recording
        [ft_sbj_newfts_ch_start_rec, stp_ftepcs_new_start_rec] = process_epochs(num_epcs_total_start_rec, ...
                                                                              rm_data, eeg_ck_start, epc_len_data, ...
                                                                              start_recording_old_fts, ch, Fs, b, a, python_directory, ...
                                                                              d_delta, d_theta, d_alpha, d_beta, d_gamma, ...
                                                                              iflinux);


        % Check if the newly generated starting points and pre-generated
        % starting points are the same (test)
        %if ~(isequal(pre_sleep_start_points, stp_ftepcs_new_pre_sleep)) 
         %   error('Precalculated starting points for unichannel features of pre-sleep period seem to differ from the one obtained during epoching for feature extraction. CHECK THIS!')
        %end 

        if ~(isequal(start_recording_start_points, stp_ftepcs_new_start_rec))
            error('Precalculated starting points for unichannel features of 5 minutes at the start of the recording seem to differ from the one obtained during epoching for feature extraction. CHECK THIS!')
        end 



        %toc
        % After parfor loop, restructure the data into multichannels way
        %mpiprofile viewer
        for i = 1:num_epcs_total_pre_sleep
            ft_sbj_newfts_pre_sleep{i}.ft_ch(ch,:) = ft_sbj_newfts_ch_pre_sleep{i};
        end

        for j = 1:num_epcs_total_start_rec
            ft_sbj_newfts_start_rec{j}.ft_ch(ch,:) = ft_sbj_newfts_ch_start_rec{j};
        end
    end

    %extract high_order features for pre-sleep period
    [fts_ho_pre_sleep, stp_ftepcs_new_ho_pre_sleep] = extract_ho_fts_from_mat(eeg_mat_pre_sleep, num_epcs_total_pre_sleep, epc_len_data, rm_data);

    %extract high_order features for 5 minutes at the start of the recording 
    [fts_ho_start, stp_ftepcs_new_ho_start] = extract_ho_fts_from_mat(eeg_mat_start_rec, num_epcs_total_start_rec, epc_len_data, rm_data);
    
    % double-check if generating points of the starting epochs are the same
    % as pre-calculated ones 
    %if ~(isequal(pre_sleep_start_points, stp_ftepcs_new_ho_pre_sleep)) 
    %        error('Precalculated starting points for high-order features of pre-sleep period seem to differ from the one obtained during epoching for feature extraction. CHECK THIS!')
    %end 

    if ~(isequal(start_recording_start_points,stp_ftepcs_new_ho_start)) 
            error('Precalculated starting points for high-order features of 5 minutes at the start of the recording seem to differ from the one obtained during epoching for feature extraction. CHECK THIS!')
    end 

    
    %stp_ftepcs_all_new{sbj_id} = stp_ftepcs_new;

    sbj_data_new.new_ft_mat_pre_sleep = ft_sbj_newfts_pre_sleep;
    sbj_data_new.ho_ft_mat_pre_sleep = fts_ho_pre_sleep;
    sbj_data_new.old_fts_pre_sleep = sbj_data.pre_sleep_old_fts;
    sbj_data_new.new_ft_mat_start = ft_sbj_newfts_start_rec;
    sbj_data_new.old_fts_start = sbj_data.start_recording_old_fts;
    sbj_data_new.ho_ft_mat_start = fts_ho_start;
    sbj_data_new.start_recording_start_points = sbj_data.start_recording_start_points;
    sbj_data_new.pre_sleep_start_points = sbj_data.pre_sleep_start_points;
    %sbj_data_new.num_ftepcs_all = num_ftepcs_all_new;
    %sbj_data_new.stp_ftepc_new = stp_ftepcs_new;
    sbj_data_new.ft_dscrp_new = ft_dscrp_new_with_ho;
    sbj_data_new.epc_to_asleep = sbj_data.epc_to_asleep;
    

     % Save the intermediate results for this subject in a separate file
    save(['new_fts_subject_', num2str(sbjs_ids(sbj_id)), '.mat'], 'sbj_data_new');
    
    disp(['Subject No.',num2str(sbjs_ids(sbj_id)),' successfully calculated'])

    
    % Display the elapsed time
    disp(['Total elapsed time for ', num2str(sbjs_ids(sbj_id)), 'subject']);
    
    % Stop the timer and get the elapsed time
    toc


end