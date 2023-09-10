function patientObj = identify_sleep_onset(patientObj,num_cont_epc_list,stage_onset, varargin)

% This function is used to identify from the given sleep scores: the
% baseline awake period, and relevant sleep onset periods (continuous stage
% N1 or N2 periods as onset, length can be self-defined)
%
% Author: Anastasia Ilina (Based on the code developed by Junheng Li)
%

%num_cont_epc - the number of continious epochs we require for the awake
%period to have?


%% Function input:
% patientObj: the patient object for processing
% num_cont_epc: array containing the number of consequtive epochs
%               required for the transition to be considered sleep onset.
%               If just 1 value is given, will record just 1 sleep onset
%               point according to the number of continious epochs defined
%               If 2 or more values -> will record a table of all different
%               sleep onsets specified by num_cont_epcs arrays

%This function basically extracts an array of indices for the first
%continious periods of chosen duration (defined by the num_cont_epcs) of
%the patient being awake and the chosen sleep stage (N1/N2), both for
%artefacts taken and not taken into account

%% Logs of code

% 15 June 2021, code started and modifying
% 21 October 2021, code modified, will run twice and return two staging
% results before and after artefact marker.

%% Varargin




if isempty(varargin)
    if_no_artefact_rejection = 0;
else
    for i = 1:length(varargin)
        switch varargin{i}
            case 'No Artefact Rejection'
                if_no_artefact_rejection = 1;
        end
    end
end


%% Discard checking
if patientObj.Sbj_discard == 1
    disp('Subject discarded, not preceding')
    return
end

%if isempty(varargin)
%    sleep_onset_type = 'Minute';


%else
%   for i = 1:length(varargin)
%        switch varargin{i}
%            case 'Python Directory'
%                python_directory = varargin{i+1};
%        end
%    end
%end


% Initialization
patientObj.if_valid_stage = 1;

%scoring = patientObj.Scoring_clean; %put the clean epoch scorings (epochs with artefacts removed from array) into a separate variable
%art_idx = patientObj.Artefact_epcs;              % Original length epochs indicating artefact locations
%     scoring = patientObj.Scoring_labels;
%     art_idx = zeros(length(scoring),1);        % Ignore artefact rejection results for now

%max_num_epc = length(scoring);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find continuous awake baseline starts


score_len = patientObj.Scoring_period;

for num_cont_epc = num_cont_epc_list

    if  if_no_artefact_rejection == 0

        clear scoring art_idx base_start Stage_start

        scoring = patientObj.Scoring_clean; %put the clean epoch scorings (epochs with artefacts removed from array) into a separate variable
        art_idx = patientObj.Artefact_epcs;              % Original length epochs indicating artefact locations
        max_num_epc = length(scoring);


        base_start = find(scoring == 0, 1); %finding the index of the first score of 0 (awake)
        if isempty(base_start)    % No such stage
            patientObj.if_valid_stage = 0;% patient is not valid for sleep scoring
            patientObj.stage_msg = 'No awake stage found';
            return
        end
        id_Onset = 1;
        max_loop = 30;
        if base_start >= max_num_epc-num_cont_epc+2 %if the first index of the awake comes before the maximum numher of scoring epochs minus the number of continipus epochs in a row that we require (i.e.2 awake epochs in a row), we discard it
            patientObj.stage_msg = 'No awake stage found';
            patientObj.if_valid_stage = 0;
            return
        end
        %do not fully understand this
        while sum(scoring(base_start:base_start+num_cont_epc-1)) ~= 0 && id_Onset<max_loop %while the sum of scorings from basline awake period to the end of the continious epochwe requre is not equal to zero and while id_Onset is samller than maxloop
            id_Onset = id_Onset+1;
            base_start = find(scoring == 0, id_Onset);%we asign new index awake epoch scorign to the basestart
            %if there are a long period of the subject being awake (stage 0), we
            %iterate through epochs to find the one where the person falls asleep
            %(stage 1) until we reach the maximal number of loops (30 iterations)
            if length(base_start)<id_Onset
                patientObj.if_valid_stage = 0;
                patientObj.stage_msg = 'No continuous awake period found';
                return
            end
            base_start = base_start(id_Onset);
            if base_start >= max_num_epc-num_cont_epc+2
                patientObj.if_valid_stage = 0;
                patientObj.stage_msg = 'No continuous awake period found';
                return
            end
        end

        idx_noart = find(art_idx == 0);    % Epochs where there are no artefacts
        patientObj.base_epcs = idx_noart(base_start:base_start+num_cont_epc-1); %we assings the baseline awaks epochs with no artefact to the patient
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the continous stage of interest start

        % Make sure that the N2 period taken are behind the awake period taken
        if base_start+num_cont_epc >= max_num_epc-num_cont_epc+2
            patientObj.if_valid_stage = 0;
            patientObj.stage_msg = 'Awake found at the end of recording';
            return
        end
        scoring = scoring(base_start+num_cont_epc:end);
        %read out input to the function (which stage we're interested in)
        if stage_onset == 'N1'
            stage_num = 1;
            sum_cont = stage_num*num_cont_epc; % as we use a sum of teh scoring as a way to check for the presence of continipuous sleeping stages, the sume of stage 2 sleep socire would be 2 times greater than for the stage 1, hence for this parameter we need to multiply it by the sleeping score values
        elseif stage_onset == 'N2'
            stage_num = 2;
            sum_cont = stage_num*num_cont_epc;
        end

        Stage_start = find(scoring == stage_num, 1); %assign the first index of the scoring array whihc contains the numbe of desired stage
        if isempty(Stage_start)    % No such stage
            patientObj.if_valid_stage = 0;
            patientObj.stage_msg = 'No sleep stage after continuous awake';
            return
        end

        if Stage_start >= max_num_epc-num_cont_epc+2
            patientObj.if_valid_stage = 0;
            patientObj.stage_msg = 'No sleep stage after continuous awake';
            return
        end

        % Find first continuous stage of num_cont_epc duration


        id_Onset = 1;
        while sum(scoring(Stage_start:Stage_start+num_cont_epc-1)) ~= sum_cont && id_Onset<max_loop %while the sum of the scores of the required numbe rof sntinious epochs is not equa to teh sumcount (parameter we defined above), and the numebr of iterations is less than teh maximal number of iterations that we defined
            id_Onset = id_Onset+1;
            Stage_start = find(scoring == stage_num, id_Onset);
            if length(Stage_start)<id_Onset
                patientObj.if_valid_stage = 0;
                patientObj.stage_msg = 'No continuous sleep stage found';
                return
            end
            Stage_start = Stage_start(id_Onset);
            if Stage_start >= max_num_epc-num_cont_epc+2
                patientObj.if_valid_stage = 0;
                patientObj.stage_msg = 'No continuous sleep stage found';
                return
            end
        end
        if sum(scoring(Stage_start:Stage_start+num_cont_epc-1)) ~= sum_cont
            patientObj.if_valid_stage = 0;
            patientObj.stage_msg = 'No continuous sleep stage found';
            return
        end


        % Correct the right index
        Stage_start = Stage_start+base_start+num_cont_epc-1;

        % Extract information about how long did the sleep initiated lasted (N1
        % included and excluded

        i= Stage_start+1;
        sleep_N1included_lenght = 1;
        sleep_N1excluded_lenght = 1;

        while scoring(i) > 0
            sleep_N1included_lenght = sleep_N1included_lenght + 1;
            if i == length(scoring)
                break
            else

                i = i +1;
            end
        end

        while i<length(scoring) || scoring(i) > 1
            sleep_N1excluded_lenght = sleep_N1excluded_lenght +1;
            if i == length(scoring)
                break
            else
                i = i +1;
            end
        end

        % Store the sleep onset(s) identified in a table that describes the
        % definition used in sleep onset identificationa and the epochs of the
        % sleep onset

        if find(num_cont_epc_list == num_cont_epc, 1) == 1 % if it's teh first onset calcualted, initialise a table
            patientObj.sleep_onset = table;
            patientObj.sleep_onset.OnsetType = strcat('First ', stage_onset, 'stage lasting for ', string(num_cont_epc*score_len/60), 'min');
            patientObj.sleep_onset.OnsetDuration = num_cont_epc;

            patientObj.sleep_onset.OnsetEpochs = {idx_noart(Stage_start:Stage_start+num_cont_epc-1)}; %exclude the artefacts
            patientObj.sleep_onset.ifCleanOnset = id_Onset == 1;
            patientObj.sleep_onset.EpochsOfSleepWithN1 = sleep_N1excluded_lenght;
            patientObj.sleep_onset.EpochsOfSleepWithoutN1 = sleep_N1excluded_lenght;
        else
            new_onset = table; % for consecutive onsets calculated, concatenate the table to the pre-existing one
            new_onset.OnsetType = strcat('First ', stage_onset, 'stage lasting for ', string(num_cont_epc*score_len/60), 'min');
            new_onset.OnsetDuration = num_cont_epc;
            new_onset.OnsetEpochs = NaN;
            new_onset.OnsetEpochs = {idx_noart(Stage_start:Stage_start+num_cont_epc-1)};
            new_onset.ifCleanOnset = id_Onset == 1;
            new_onset.EpochsOfSleepWithN1 = sleep_N1excluded_lenght;
            new_onset.EpochsOfSleepWithoutN1 = sleep_N1excluded_lenght;
            patientObj.sleep_onset = [patientObj.sleep_onset; new_onset];


        end

        patientObj.stage_msg = 'All good';

    end

    %% Running again for scorings without artefact markers

    clear scoring art_idx base_start Stage_start

    % Initialization
    patientObj.if_valid_stage_noart = 1;

    scoring = patientObj.Scoring_labels; %let's take out the sleep staging information (720 scores for each 30 s epoch) into a separate variable
    fake_art_idx = zeros(length(scoring),1);        % Ignore artefact rejection results for now

    max_num_epc = length(scoring);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find continuous awake baseline starts
    base_start = find(scoring == 0, 1); %return the first index of the sleep scoring array, which has 0 (awake stage)
    if isempty(base_start)    % If no awake stage has been found (the base_start is empty)
        patientObj.if_valid_stage_noart = 0; %the object is not valid for sleep stagin
        patientObj.stage_msg_noart = 'No awake stage found';
        return
    end
    id_Onset = 1;
    max_loop = 30;
    if base_start >= max_num_epc-num_cont_epc+2
        patientObj.stage_msg_noart = 'No awake stage found';
        patientObj.if_valid_stage_noart = 0;
        return
    end
    while sum(scoring(base_start:base_start+num_cont_epc-1)) ~= 0 && id_Onset<max_loop
        id_Onset = id_Onset+1;
        base_start = find(scoring == 0, id_Onset);
        if length(base_start)<id_Onset
            patientObj.if_valid_stage_noart = 0;
            patientObj.stage_msg_noart = 'No continuous awake period found';
            return
        end
        base_start = base_start(id_Onset);
        if base_start >= max_num_epc-num_cont_epc+2
            patientObj.if_valid_stage_noart = 0;
            patientObj.stage_msg_noart = 'No continuous awake period found';
            return
        end
    end

    idx_no_art_rejection = find(fake_art_idx == 0);    % Epochs where there are no artefacts
    patientObj.base_epcs_noart = idx_no_art_rejection(base_start:base_start+num_cont_epc-1);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the continous stage start

    % Make sure that the N2 period taken are behind the awake period taken
    if base_start+num_cont_epc >= max_num_epc-num_cont_epc+2
        patientObj.if_valid_stage_noart = 0;
        patientObj.stage_msg_noart = 'Awake found at the end of recording';
        return
    end
    scoring = scoring(base_start+num_cont_epc:end);

    if stage_onset == 'N1'
        stage_num = 1;
        sum_cont = stage_num*num_cont_epc;
    elseif stage_onset == 'N2'
        stage_num = 2;
        sum_cont = stage_num*num_cont_epc;
    end

    Stage_start = find(scoring == stage_num, 1);
    if isempty(Stage_start)    % No such stage
        patientObj.if_valid_stage_noart = 0;
        patientObj.stage_msg_noart = 'No sleep stage after continuous awake';
        return
    end
    id_Onset = 1;
    max_loop = 30;
    if Stage_start >= max_num_epc-num_cont_epc+2
        patientObj.if_valid_stage_noart = 0;
        patientObj.stage_msg_noart = 'No sleep stage after continuous awake';
        return
    end

    % Find first continuous stage minute
    while sum(scoring(Stage_start:Stage_start+num_cont_epc-1)) ~= sum_cont && id_Onset<max_loop
        id_Onset = id_Onset+1;
        Stage_start = find(scoring == stage_num, id_Onset);
        if length(Stage_start)<id_Onset

            patientObj.if_valid_stage_noart = 0;
            patientObj.stage_msg_noart = 'No continuous sleep stage found';
            return
        end
        Stage_start = Stage_start(id_Onset);
        if Stage_start >= max_num_epc-num_cont_epc+2
            patientObj.if_valid_stage_noart = 0;
            patientObj.stage_msg_noart = 'No continuous sleep stage found';
            return
        end
    end
    if sum(scoring(Stage_start:Stage_start+num_cont_epc-1)) ~= sum_cont
        patientObj.if_valid_stage_noart = 0;
        patientObj.stage_msg_noart = 'No continuous sleep stage found';
        return
    end


    % Correct the right index
    Stage_start = Stage_start+base_start+num_cont_epc-1;

    % Extract information about how long did the sleep initiated lasted (N1
    % included and excluded

    i= Stage_start+1;
    sleep_N1included_lenght = 1;
    sleep_N1excluded_lenght = 1;
    while scoring(i) > 0
        sleep_N1included_lenght = sleep_N1included_lenght +1;
        i = i +1;
    end

    while scoring(i) > 1
        sleep_N1excluded_lenght = sleep_N1excluded_lenght +1;
        i = i +1;
    end


    % Store the sleep onset(s) identified in a table that describes the
    % definition used in sleep onset identificationa and the epochs of the
    % sleep onset

    if find(num_cont_epc_list == num_cont_epc, 1) == 1 % if it's teh first onset calcualted, initialise a table
        patientObj.sleep_onset_noart = table;
        patientObj.sleep_onset_noart.OnsetType = strcat('First ', stage_onset, 'stage lasting for ', string(num_cont_epc*score_len/60), 'min');
        patientObj.sleep_onset_noart.OnsetDuration = num_cont_epc;
        patientObj.sleep_onset_noart.OnsetEpochs = {idx_no_art_rejection(Stage_start:Stage_start+num_cont_epc-1)}; %exclude the artefacts
        patientObj.sleep_onset_noart.ifCleanOnset = id_Onset == 1;
        patientObj.sleep_onset_noart.EpochsOfSleepWithN1 = sleep_N1excluded_lenght;
        patientObj.sleep_onset_noart.EpochsOfSleepWithoutN1 = sleep_N1excluded_lenght;
    else
        new_onset_noart = table; % for consecutive onsets calculated, concatenate the table to the pre-existing one
        new_onset_noart.OnsetType = strcat('First ', stage_onset, 'stage lasting for ', string(num_cont_epc*score_len/60), 'min');
        new_onset_noart.OnsetDuration = num_cont_epc;
        new_onset_noart.OnsetEpochs = {idx_no_art_rejection(Stage_start:Stage_start+num_cont_epc-1)};
        new_onset_noart.ifCleanOnset = id_Onset == 1;
        new_onset_noart.EpochsOfSleepWithN1 = sleep_N1excluded_lenght;
        new_onset_noart.EpochsOfSleepWithoutN1 = sleep_N1excluded_lenght;
        patientObj.sleep_onset_noart = [patientObj.sleep_onset_noart; new_onset_noart];

    end

    patientObj.stage_msg_noart = 'All good';
end