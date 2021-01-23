% Run this file after Convert2Scherbaum
cd 'C:\Users\dpen466\Google Drive\Phd (1)\Experiment\Ecological_psych\Erb2016_analysis/AnalysisN149/Spilt/Response_Repeat/RR_all'
% below identify all trials with same trial values in x, y, t
% then concatenate them into cell arrays and then into the same format as
% Scherbaum's data structure with the desired predictors and save it in
% their separate files for each participant in response repeat and response switch. 
fieldnames= {'congruency'; 'previous_congruency'; 'response'; 'rt'; 'x'; 'y'; 't'}; % change if predictor changes
fieldnames2= {'congruency'; 'response'; 'rt'; 'x'; 'y'; 't'};
for participant_n = 1: length(file_names)
    trial_uniq = unique(dynamics(participant_n).Trial);% Find the unique trials number in dyanmics file
    x = cell(length(allTrials(participant_n).Trial), 1);% generate cell arrays of 0 to be filled in
    y = cell(length(allTrials(participant_n).Trial), 1);
    t = cell(length(allTrials(participant_n).Trial), 1);
    fprintf('Spilting participant: %d\n', participant_n)
    for trial_n = 1:length(trial_uniq)
        x_data = dynamics(participant_n).x(dynamics(participant_n).Trial==trial_uniq(trial_n));%code to extract x, y and t from the correct trial
        y_data = dynamics(participant_n).y(dynamics(participant_n).Trial==trial_uniq(trial_n));
        t_data = dynamics(participant_n).t(dynamics(participant_n).Trial==trial_uniq(trial_n));
        x{trial_n,1} = x_data';% start filling in with double arrays
        y{trial_n,1} = y_data';
        t{trial_n,1} = t_data';
        congruency= num2cell(allTrials(participant_n).Congruency); % convert predictors from double to cell
        pre_cong =num2cell(allTrials(participant_n).Previous_congruency);
        trial_sequence=num2cell(allTrials(participant_n).trial_sequence);
        response=num2cell(allTrials(participant_n).LocationTouched);
        rt= num2cell(allTrials(participant_n).RT);
        data =[congruency pre_cong trial_sequence response rt x y t]; % binds all variables into one cell array
    end
    % this might be problematic
    data_RR = data([data{:,3}] == 1,:);
    data_RR= data_RR (:,[1,2,4:8]);
    data_RS = data([data{:,3}] == 2,:);
    data_RS= data_RS (:,[1,2,4:8]);
    data_RRPC =data_RR([data_RR{:,2}] == 1,:);
    data_RRPC= data_RRPC (:,[1,3:7]);
    % data_RRIC= data_RR([data_RR{:,1}] == 2,:); %for spilting by current
    %congruency
    %data_RRIC= data_RRIC (:,[2:7]);
    data_RSPC = data_RS([data_RS{:,2}] == 1,:);
    data_RSPC = data_RSPC (:,[1,3:7]);
    %save response repeat that include both current incongruent and
    %congruent trials for each participants in RR_all
    cd  '../../Response_Repeat/RR_all'
    trials = cell2struct(data_RR, fieldnames, 2)';
    fname = sprintf('log%d.mat', participant_n);
    save(fname, 'trials'); %save the file as logn.mat
    %do the same but only for previous congruent trials
    cd '../RR_PC' 
    trials = cell2struct(data_RRPC, fieldnames2, 2)';
    fname= sprintf('log%d.mat', participant_n);
    save(fname, 'trials'); %save the file as logn.mat
    %do the same but only for current Incongruent trials
    %cd '../RR_Incongruent' 
    %trials = cell2struct(data_RRIC, fieldnames2, 2)';
    %fname= sprintf('log%d.mat', participant_n);
    %save(fname, 'trials'); %save the file as logn.mat
    
    % save files for response switch
    cd  '../../Response_Switch/RS_all'
    trials = cell2struct(data_RS, fieldnames, 2)';
    fname= sprintf('log%d.mat', participant_n);% convert the cell array to struct in the same format as Sherbaum's data
    save(fname, 'trials'); % save the file as logn.mat
    cd '../RS_PC' 
    trials = cell2struct(data_RSPC, fieldnames2, 2)';
    fname= sprintf('log%d.mat', participant_n);
    save(fname, 'trials'); % save the file as logn.mat
end
clear

