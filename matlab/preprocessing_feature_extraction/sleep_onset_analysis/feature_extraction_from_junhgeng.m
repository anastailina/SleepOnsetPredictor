%% Preparations

% Subject group
cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All'
load fteval_EWSallsbj_20Dec.mat
sbj_id = 1;

sbjs_ids = sbj_remained;



%% Continuous feature evaluation pipelines

% Original analysis parameters
Stage = 'N2';    % Stage onset to study
num_sbjs = length(sbjs_ids);

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
cd '/mnt/LongTermStorage/MESA Sleep-Onset Extracted All/anastasia_analysis'
for sbj_id = 1:num_sbjs
    
    % Start the timer
    tic;

    % Feature evaluation
    % load the eeg data for the subject
    eeg_mat_now = eeg_asleepper_all{sbj_id}.ch;
    num_chs = max(size(eeg_mat_now));
    
    % load the starting points of the epoch from Junheng data
    stp_points_sbj = stp_ftepcs_all{sbj_id};
    
    % load the features extracted by Junheng 
    ft_one_sbj = ft_sbj{sbj_id};

    % initialise counter for the channel index
    idx_realch = 0;

    % initialise cell to store fts for extracted samples 
    ft_sbj_newfts = cell(1, num_epcs_total);

    % iterate through the channel
    for ch = 1:num_chs

        % increment the channel counter
        idx_realch = idx_realch+1;
        
        % get the channel of our the eeg matrix
        eeg_ck = eeg_mat_now{idx_realch};
        
        len_now = length(eeg_ck);
        epc_len_data = epc_len*Fs;
        ovlap_data = floor(epc_len_data*ovlap);
        rm_data = epc_len_data-ovlap_data;
        
        % calculate the total number of epochs given overlapping
        num_epcs_total = floor((len_now - ovlap_data)/rm_data);  

        % increment the total number fo the epochs that the subject has
        num_ftepcs_all_new(sbj_id) = num_ftepcs_all_new(sbj_id) + num_epcs_total; 
        
        %initialise the variable storing new starting points storage for
        %each of the epochs
        stp_ftepcs_new = zeros(num_epcs_total,1);

        % Initialize the output cell array outside the parfor loop
        ft_sbj_newfts_ch = cell(1, num_epcs_total);
        
        % iterate through all of the epochs 
        for jj = 1:num_epcs_total
            % Starting timepoint of current epoch
            stp = (jj-1) * rm_data + 1; 
            stp_ftepcs_new(jj) = stp;

            % Current data chunk
            data_unit = eeg_ck(stp:stp+epc_len_data-1);  
            
            % 0.5 is technical treshold of EEG. If subject is moving
            % significantly, it reaches 0.5 value -> movement artefact ->
            % discard the epoch
            if max(data_unit) == 0.5 
                ft_sbj_newfts_ch{jj} = NaN(1, 47);
            
            else
                % load the samples from junheng

                %get the index of the sampling (he sampled with 0.5
                %overlap, i'm doing 0 overlap)
                idx_ft_sample = 2 * (jj - 1) + 1;

                %get the old features for the channel
                fts_old = ft_one_sbj.ck{idx_ft_sample}.ft_ch(ch,:);
                
                % get delta, theta, alpha and beta powers
                delta_power = fts_old(1);
                theta_power = fts_old(2);
                alpha_power = fts_old(3);
                beta_power = fts_old(4);
                
                % extract new features
                [ftt_new, ft_dscrp_new] = sleep_ft_extract_optimised(data_unit, ...
                    Fs, delta_power, theta_power, alpha_power, beta_power, ...
                    b, a, python_directory, d_delta, d_theta, d_alpha, d_beta, ...
                    d_gamma, iflinux);
                
                % convert the structure to cell
                ftt1_new = struct2cell(ftt_new);

                % convert the cell to mat and store it in the cell allong
                % with all other channel samples
                ft_sbj_newfts_ch{jj} = transpose(cell2mat(ftt1_new));
            end
        end
        
        % After parfor loop, restructure the data into multichannels way
        
        for jj = 1:num_epcs_total
            ft_sbj_newfts{jj}.ft_ch(ch,:) = ft_sbj_newfts_ch{jj};
        end
    end
    
    %stp_ftepcs_all_new{sbj_id} = stp_ftepcs_new;

    sbj_data_new = 

     % Save the intermediate results for this subject in a separate file
    save(['new_fts_subject_', num2str(sbj_id), '_results.mat'], 'ft_sbj_newfts', 'num_ftepcs_all_new', 'stp_ftepcs_new');
    
    disp(['Subject No.',num2str(sbj_id),' successfully calculated'])

    % Stop the timer and get the elapsed time
    elapsedTime = toc;

    % Display the elapsed time
    disp(['Total elapsed time: ', num2str(elapsedTime), ' seconds']);


end