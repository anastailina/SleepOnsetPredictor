function [ft_sbj_newfts_ch, stp_ftepcs_new] = process_epochs(num_epcs_total, rm_data, eeg_ck, epc_len_data, ft_one_sbj, ch, Fs, b, a, python_directory, d_delta, d_theta, d_alpha, d_beta, d_gamma, iflinux)

   % % Set up Python
    %python_directory = '/home/anastasia/Python-3.10.11/python';
    %iflinux = 1;
    %setup_python(python_directory, iflinux);
% 

    % initialise the variable storing new starting points storage for each of the epochs
    stp_ftepcs_new = zeros(num_epcs_total,1);

    % Initialize the output cell array outside the parfor loop
    ft_sbj_newfts_ch = cell(1, num_epcs_total);
     
    %mpiprofile on
    %tic;

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
            ft_sbj_newfts_ch{jj} = NaN(1, 19);
        else
            % get the index of the sampling (he sampled with 0.5
            % overlap, i'm doing 0 overlap)
           % idx_ft_sample = 2 * (jj - 1) + 1;

           % for already downsampled
            idx_ft_sample = jj; 

            % get the old features for the channel
            fts_old = ft_one_sbj{idx_ft_sample}.ft_ch(ch,:);

            % get delta, theta, alpha and beta powers
            delta_power = fts_old(1);
            theta_power = fts_old(2);
            alpha_power = fts_old(3);
            beta_power = fts_old(4);

            % extract new features
            ftt_new = sleep_ft_extract_optimised(data_unit, ...
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
    %mpiprofile viewer
    %toc
end