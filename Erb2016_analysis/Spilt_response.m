%%%
% README: (Optional) To be run after Convert2Scherbaum.m. This script spilt responses according to the analyst's choosing 
% and convert to the appropriate data format for TCMR. In this particular
% script we are spilting the data according to previous_congruency so we
% can run TCMR separately for each trial sequence type (i.e., response
% repeat
%%%

%% Processing response repeat and response switch
% Navigate to folder where Response repeat data file will be stored
cd 'C:\Users\dpen466\Google Drive\Phd (1)\Share_with_Chris\TCMR_Analysis_Pipeline\Erb2016_analysis\AnalysisN20\Spilt\Response_Repeat\RR_all' 
% Below identify all trials with same trial values in x, y, t
% then concatenate them into cell arrays and then into the same format as
% Scherbaum's data structure with the desired predictors and save it in
% their separate files for each participant in response repeat and response switch. 
fieldnames= {'congruency'; 'previous_congruency'; 'response'; 'rt'; 'x'; 'y'; 't'}; % change if predictor changes
for participant_n = 1: length(file_names)
    trial_uniq = unique(dynamics(participant_n).Trial); % Find the unique trials number in dyanmics file
    x = cell(length(allTrials(participant_n).Trial), 1); % generate cell arrays of 0 to be filled in
    y = cell(length(allTrials(participant_n).Trial), 1);
    t = cell(length(allTrials(participant_n).Trial), 1);
    fprintf('Spilting participant: %d\n', participant_n)
    for trial_n = 1:length(trial_uniq)
        x_data = dynamics(participant_n).x(dynamics(participant_n).Trial==trial_uniq(trial_n));% code to extract x, y and t from the correct trial
        y_data = dynamics(participant_n).y(dynamics(participant_n).Trial==trial_uniq(trial_n));
        t_data = dynamics(participant_n).t(dynamics(participant_n).Trial==trial_uniq(trial_n));
        x{trial_n, 1} = x_data';% start filling in with double arrays
        y{trial_n, 1} = y_data';
        t{trial_n, 1} = t_data';
        congruency = num2cell(allTrials(participant_n).Congruency); % convert predictors from double to cell
        pre_cong = num2cell(allTrials(participant_n).Previous_congruency);
        trial_sequence = num2cell(allTrials(participant_n).trial_sequence);
        response = num2cell(allTrials(participant_n).LocationTouched);
        rt = num2cell(allTrials(participant_n).RT);
        data = [congruency pre_cong trial_sequence response rt x y t]; % binds all variables into one cell array
    end
    % Spilt data file accoridng to trial sequence
    data_RR = data([data{:,3}] == 1,:);% Response repeat
    data_RR = data_RR (:,[1,2,4:8]);
    data_RS = data([data{:,3}] == 2,:);% Response swicth
    data_RS = data_RS (:,[1,2,4:8]);
    
    % Save response spilt data file to their respective folders
    % Response repeat
    cd  '../../Response_Repeat/RR_all'
    trials = cell2struct(data_RR, fieldnames, 2)';
    fname = sprintf('log%d.mat', participant_n);
    save(fname, 'trials'); 
    % Reponse switch
    cd  '../../Response_Switch/RS_all'
    trials = cell2struct(data_RS, fieldnames, 2)';
    fname= sprintf('log%d.mat', participant_n);
    save(fname, 'trials');
end
%%
cd ../../../../
clear

