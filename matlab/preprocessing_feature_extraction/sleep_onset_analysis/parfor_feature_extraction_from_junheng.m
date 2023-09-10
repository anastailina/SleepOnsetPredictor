%% Preparations

% Subject group
clear all
cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All/anastasia_analysis'

load sbj_remained_30mins.mat
load ft_dscrp_new.mat


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

%for sbj_id = 1:num_sbjs

    % Start the timer
    tic;
    sbj_id = 1;
    % load the subject
    load(['subject_', num2str(sbjs_ids(sbj_id)), '.mat'])
    

    % Feature evaluation
    % load the eeg data for the subject
    eeg_mat_now = sbj_data.eeg_mat;
    num_chs = max(size(eeg_mat_now));
    
    % load the starting points of the epoch from Junheng data
    stp_points_sbj = sbj_data.stp_points;
    
    % load the features extracted by Junheng 
    ft_one_sbj = sbj_data.ft_mat.ck;

    % initialise counter for the channel index
    idx_realch = 0;

    % get the channel of our the eeg matrix    
       
    len_now = length(eeg_mat_now{1});
    epc_len_data = epc_len*Fs;
    ovlap_data = floor(epc_len_data*ovlap);
    rm_data = epc_len_data-ovlap_data;
        
    % calculate the total number of epochs given overlapping
    num_epcs_total = floor((len_now - ovlap_data)/rm_data);  

    % initialise cell to store fts for extracted samples 
    ft_sbj_newfts = cell(1, num_epcs_total);

    % iterate through the channel
    for ch = 1:num_chs

        % increment the channel counter
        idx_realch = idx_realch+1;
        
        % get the channel of our the eeg matrix
        eeg_ck = eeg_mat_now{idx_realch};

         % increment the total number of the epochs that the subject has
        num_ftepcs_all_new(sbj_id) = num_ftepcs_all_new(sbj_id) + num_epcs_total; 


        [ft_sbj_newfts_ch, stp_ftepcs_new, num_ftepcs_all_new] = process_epochs(sbj_id, num_epcs_total, rm_data, eeg_ck, epc_len_data, ft_one_sbj, ch, Fs, b, a, python_directory, d_delta, d_theta, d_alpha, d_beta, d_gamma, iflinux, num_ftepcs_all_new);
        %toc
        % After parfor loop, restructure the data into multichannels way
        %mpiprofile viewer
        for jj = 1:num_epcs_total
            ft_sbj_newfts{jj}.ft_ch(ch,:) = ft_sbj_newfts_ch{jj};
        end
    end
    
    %stp_ftepcs_all_new{sbj_id} = stp_ftepcs_new;

    sbj_data_new.ft_mat_new = ft_sbj_newfts;
    sbj_data_new.num_ftepcs_all = num_ftepcs_all_new;
    sbj_data_new.stp_ftepc_new = stp_ftepcs_new;
    sbj_data_new.ft_dscrp_new = ft_dscrp_new;
    sbj_data_new.epc_to_asleep = sbj_data.epc_to_asleep;
    

     % Save the intermediate results for this subject in a separate file
    save(['new_for_fts_subject_', num2str(sbjs_ids(sbj_id)), '.mat'], 'sbj_data_new');
    
    disp(['Subject No.',num2str(sbjs_ids(sbj_id)),' successfully calculated'])

    
    % Display the elapsed time
    disp(['Total elapsed time for ', num2str(sbjs_ids(sbj_id)), 'subject']);
    
    % Stop the timer and get the elapsed time
    toc


%end