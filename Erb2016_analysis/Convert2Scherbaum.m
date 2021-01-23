



% load csv files from cleaned files of dynamics and Alltrials
addpath(genpath('Source'))
cd 'C:\Users\dpen466\Google Drive\Phd (1)\Share_with_Chris\TCMR_Analysis_Pipeline\Erb2016_data\All_trials_cleaneddata'
clear
% All_trials file
file_names = dir ('*.csv');
file_names = struct2cell(file_names);
file_names = natsortfiles(file_names(1,:));%code to sort the loaded file name into the right order
% here I am converting all the data from dynamics and Alltrials into  their
% respective struct with each row representing one participant
for i = 1:length(file_names)
    disp(file_names(i))
    allTrials(i)=csv2struct(char(file_names(i)));
end
% Dynamics file
cd '../dynamics_cleaneddata'
file_names = dir ('*.csv');
file_names = struct2cell(file_names);
file_names = natsortfiles(file_names(1,:));
for i = 1:length(file_names)
    disp(file_names(i))
    dynamics(i) = csv2struct(char(file_names(i)));
end

cd '../../Erb2016_analysis/AnalysisN149/Include_RR_RS/'
% below identify all trials with same trial values in x, y, t
% then concatenate them into cell arrays and then into the same format as
% Scherbaum's data structure with the desired predictors and save it in
% their separate files for each participant. 
fieldnames= {'congruency'; 'response'; 'previous_congruency';'rt';'x';'y';'t'}; % change if predictor changes
for participant_n = 1:length(file_names) 
    trial_uniq = unique(dynamics(participant_n).Trial);% Find the unique trials number in dyanmics file
    x = cell(length(allTrials(participant_n).Trial), 1);% generate cell arrays of 0 to be filled in
    y = cell(length(allTrials(participant_n).Trial), 1);
    t = cell(length(allTrials(participant_n).Trial), 1);
    fprintf('Processing participant: %d\n', participant_n)
    for trial_n = 1:length(trial_uniq)
        x_data = dynamics(participant_n).x(dynamics(participant_n).Trial == trial_uniq(trial_n));% code to extract x, y and t from the correct trial
        y_data = dynamics(participant_n).y(dynamics(participant_n).Trial == trial_uniq(trial_n));
        t_data = dynamics(participant_n).t(dynamics(participant_n).Trial == trial_uniq(trial_n));
        x{trial_n, 1} = x_data';% start filling in with double arrays
        y{trial_n, 1} = y_data';
        t{trial_n, 1} = t_data';
        congruency = num2cell(allTrials(participant_n).Congruency); % convert predictors from double to cell
        pre_cong =num2cell(allTrials(participant_n).Previous_congruency);
        response = num2cell(allTrials(participant_n).LocationTouched);
        rt = num2cell(allTrials(participant_n).RT);
        data = [congruency response pre_cong rt x y t]; % binds all variables into one cell array
    end
    trials = cell2struct(data, fieldnames, 2)'; % convert the cell array to struct in the same format as Sherbaum's data
    fname = sprintf('log%d.mat', participant_n);
    save(fname, 'trials'); % save the file as logn.mat
end
cd ../
% right now the logn file saved, the n number for each participant does not
% match the once in raw data, so removed participants are skip, i.e. if
% participant 4 is removed then log4 is actually participant 4, perhaps need
% to write a code so that logn reflect the actual participant number?
%Spilt_response

