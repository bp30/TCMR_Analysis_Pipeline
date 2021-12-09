
%%%
% README: Script to convert cleaned csv files produced by CleanfilesFlanker.R to
% the appropriate data format for TCMR. It take the cleaned All_trials and
% dynamics cvs files to logn.m files that can then be analysed with
% Scherbaum's codes
%%%

%% Add path to scripts required
addpath(genpath('Source'))
clear

%% Processing All_trials files
% Navigate to folder containing cleaned data files for All_trials 
cd '../Erb2016_data/All_trials_cleaneddata/'
% Identify all the All_trial files in the folder
file_names = dir ('*.csv');
% Convert the struct of filenames for Alltrials into cells
file_names = struct2cell(file_names);
% Extract the file names and sort then by order so now each column is one
% participant
file_names = natsortfiles(file_names(1,:));
% Using file_names than extract data from each Alltrials file and append it
% to a struct
for i = 1 : length(file_names)
    fprintf('Loading All_trials for participant: %d\n', i)
    allTrials(i) = csv2struct(char(file_names(i)));
end

%% Processing of dynamics files
% Navigate to the directory
cd '../dynamics_cleaneddata'
% Same process as All_trials csv files
file_names = dir ('*.csv');
file_names = struct2cell(file_names);
file_names = natsortfiles(file_names(1,:));
for i = 1:length(file_names)
    fprintf('Loading dynamics for participant: %d\n', i)
    dynamics(i) = csv2struct(char(file_names(i)));
end

%% Convert to data format suitable for TCMR
% Navigate to folder where the analysis script is located
cd '../../Erb2016_analysis/AnalysisN20/Include_RR_RS/'

% The codes below identify all trials with same trial number for x, y, t
% then concatenate them into cell arrays and then into the same format as
% TCMR data structure with the desired predictors and save it separately
% for each participant as logn.m

% Variable names that need to be extracted
fieldnames = {'congruency'; 'response'; 'previous_congruency';'rt';'x';'y';'t'}; % change if predictor changes
% Extract the necessary information and save to logn.m file for each
% participant
for participant_n = 1:length(file_names) 
    trial_uniq = unique(dynamics(participant_n).Trial);% Find the unique trials number in dyanmics file
    x = cell(length(allTrials(participant_n).Trial), 1);% generate cell arrays of 0 to be filled in
    y = cell(length(allTrials(participant_n).Trial), 1);
    t = cell(length(allTrials(participant_n).Trial), 1);
    fprintf('Processing participant: %d\n', participant_n)
    for trial_n = 1 : length(trial_uniq)
        % Extract x, y, and t for each trial from dynamics 
        x_data = dynamics(participant_n).x(dynamics(participant_n).Trial == trial_uniq(trial_n));
        y_data = dynamics(participant_n).y(dynamics(participant_n).Trial == trial_uniq(trial_n));
        t_data = dynamics(participant_n).t(dynamics(participant_n).Trial == trial_uniq(trial_n));
        x{trial_n, 1} = x_data';% start filling in with double arrays
        y{trial_n, 1} = y_data';
        t{trial_n, 1} = t_data';
    end   
    % Extract predictor information from All_trials
    congruency = num2cell(allTrials(participant_n).Congruency); % convert predictors from double to cell
    pre_cong = num2cell(allTrials(participant_n).Previous_congruency);
    response = num2cell(allTrials(participant_n).LocationTouched);
    rt = num2cell(allTrials(participant_n).RT);
        
    % Bind all variables into one cell array
    data = [congruency response pre_cong rt x y t]; 
    
    % Convert the cell array to struct so it is in the appropriate formate for TCMR
    trials = cell2struct(data, fieldnames, 2)'; 
    fname = sprintf('log%d.mat', participant_n);
    % Save the file as logn.mat
    save(fname, 'trials');
end
cd ../../

% Optional: To run if spilting data file is necessary
%Spilt_response


%% Note: right now the logn file saved, the n number for each participant does not
% match the once in raw data, so removed participants are skip, i.e. if
% participant 4 is removed then log4 is actually participant 5.
% Spilt_response

