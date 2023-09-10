function sampleObj_coll = ft_eval_sleep_onset(sampleObj_coll,txt_feature,ifcardiac, varargin)

% This function is used to evaluate the features and turn it into a feature
% vector with corresponding descriptions
%
% Author: Junheng Li, Anastasia Ilina

%% Log of code

% 29/06/2021, code started
% 21/10/2021, code contined to finish
% 09/05/2023, Anastasia: changed the range of the band-pass filter

%% 
if isempty(varargin)
    python_directory = '/home/anastasia/Python-3.10.11/python';

else
    for i = 1:length(varargin)
        switch varargin{i}
            case 'Python Directory'
                python_directory = varargin{i+1};
        end
    end 
end

%% Loading feature description file



fid = fopen(txt_feature,'rt');
A = textscan(fid,'%s%s');% extract the names of all of teh features we want from a text file (contains two  cells) 

ft_names = A{1}; %extract the names of the features form the first cell of A
num_ft = length(ft_names); %define the number of features we are interested in 

%% Feature evaluations

num_samps = length(sampleObj_coll);        % Total number of samples in the collection

% Enables parallel computation 
parfor i = 1:num_samps %iterating through the number of samples 

    samp_now = sampleObj_coll{i}; %extract the ith sample structure

    samp_now.ft_quality = 1; % initialise the feature quality as good, if smth is wrong -> change it later
    
    y = samp_now.data_EEG; 
    if sum(isnan(y)) == 0
        Fs = samp_now.Fs; 
        [ft,~] = sleep_ft_extract_new(y,Fs, python_directory, 1, 'Min Freq', 0.1, 'Max Freq', 100);        
        ft = struct2cell(ft); 
        ft = transpose(cell2mat(ft)); 
    else
        ft = NaN(1,num_ft);
        samp_now.ft_vec = ft;
        samp_now.ft_quality = 0;
        disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' contains NaN EEG data'])
        
    end
        
    if ifcardiac              % Possibly add ECG feature evaluation in the future
        
    end
    
    %Checking if the subject has the right snumber of the extracted
    %features 
    if length(ft) ~= num_ft
        samp_now.ft_quality = 0;
        warning(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' does not have the right feature number'])
        
    else
        %check if any extracted features of this subject contain NaN values
        samp_now.ft_vec = ft;
        if sum(isnan(ft))>0
            samp_now.ft_quality = 0;
            disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' have NaN feature value'])
            
        end

        samp_now.ft_dscrp = ft_names;
        %If all is okay,signifiy it by the ft_quality variable = 1 and
        %print out that ft-eval went well 
        %disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' has been successfully calculated'])
        
    end
    
    sampleObj_coll{i} = samp_now; %input the extracted features and corresponding iformation back into teh samples form the subject structure 

end










