function [art_all, art_ch, proportion_epochs_remaining] = rms_art_rm_mask(eeg_mat_now, epc_len, num_epc, Fs, iter, th_coeff, method)

% This function is used for general artefact rejection based on large
% epochs for removing extraordinary artefacts such as movements;
%
% Author: Anastasia Ilina, Junheng Li
%
% method: either 'RMS' or 'Otsu' for selecting the thresholding technique


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-calculation
epc_pt = epc_len*Fs;
num_eeg_ch = min(size(eeg_mat_now));
art_all = zeros(1,num_epc);
art_ch = zeros(num_eeg_ch,num_epc);

% Artefact detection
num_epcs_rejected = 0;

for ch = 1:num_eeg_ch
    
    if isempty(eeg_mat_now(ch, :))
        art_ch(ch,:) = ones(1,num_epc);                        % If this channel is discarded then no detection should be done
        continue
    end

    EEG_ch = eeg_mat_now(ch,:);   % EEG trace of current channel;
    rms_all_ch = zeros(num_epc,1);
    start_epc = zeros(num_epc,1);

    for i = 1:num_epc
        start = 1+(i-1)*epc_pt;
        eeg_epc = EEG_ch(start:start+epc_pt-1);

        rms_all_ch(i) = rms(eeg_epc);
        start_epc(i) = start;
    end
    
    % Number of iterations
    for ite = 1: iter
        if strcmp(method, 'RMS')
            med_rms = median(rms_all_ch,'omitnan');
            th = th_coeff*med_rms;
        elseif strcmp(method, 'Otsu')
        % Remove NaN values, then rescale to [0, 1]
        valid_rms = rms_all_ch(~isnan(rms_all_ch));
        scaled_rms = (valid_rms - min(valid_rms)) / (max(valid_rms) - min(valid_rms));
        otsu_level = graythresh(scaled_rms);
        th = otsu_level * (max(valid_rms) - min(valid_rms)) + min(valid_rms);
        else
            error('Unknown method. Use "RMS" or "Otsu".');
        end

        for i = 1:num_epc
            if rms_all_ch(i)> th
                num_epcs_rejected = num_epcs_rejected + 1;
                art_all(i) = 1;
                art_ch(ch,i) = 1;
                rms_all_ch(i) = NaN;
            end
        end
    end
end

proportion_epochs_remaining = 1 - sum(art_all)/num_epc;
disp(['% Epochs remaining after rejecting ones ', num2str(th_coeff), ' times the thresholding method (', method, '): ', num2str(proportion_epochs_remaining)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
