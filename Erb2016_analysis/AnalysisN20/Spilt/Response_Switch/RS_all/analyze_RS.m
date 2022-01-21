% DemoScript for TCMR Toolbox by Stefan Scherbaum (C) 2017
% Reads a subset of the data from Scherbaum et al., Cognition, 2010, study
% 1 and analyzes it using TCMR Toolbox
%
% Author: Stefan Scherbaum, University of Dresden, 2017
%
% Copyright (C) 2017 Stefan Scherbaum, stefan.scherbaum@tu-dresden.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Revision 1.0 schtefan

%% load data
% BP comments: the below codes load the data (logn.mat files; contains information on RT, congruencyN, responseN, xy coordinates and time points (t)) into the workspace and create 
% a struct (can contain numbers and strings) for each participant with basic properties of the trials (response, congruency as vectors). 
% BP comment: it generates combine each struct to a larger struct where
% each row is one participant and each column is a variablead
clear 
addpath(genpath('../../../Source'))
for p=1:20 %loop from 1 to N number of participants
    out('Loading log ',p);%BP:output on terminal saying which dataset is loading now
    load(out('./log',p,'.mat')); %BP: load logn.mat file (datafile for each participant)
    pstruct=catTrialLog(trials,{'response','congruency','previous_congruency', 'rt'});%extract basic data from each trial struct into matrices, Note rt is not need to calculate anything subsequent analysis
    pstruct.trials=trials;%add trial struct for further processing
    pdata(p)=pstruct; %bp: combine pstructs for each participant to struct file.
end

%% perform preprocessing
sf=95;%sampling was at 95 Hz
warp_samples=100;%100 time slices for time normalized trajectories;
stimlock_samples=ceil(sf*2.5);%2.5 seconds max response time, 95 samples per second --> maximum samples planned for stimulus locked trajetories - rest is filled with nan
screen_width=1280;%width of screen in the experiment

for p=1:length(pdata)
    out('Preprocessing dataset ',p,'...');
    %initialize matrices for stimulus locked raw data 
    pdata(p).x_stimlock=nan(length(pdata(p).trials),stimlock_samples); %bp:generate nan arrays for x, t and y in pdata for each participants that sets max RT to 2.5s 
    pdata(p).t_stimlock=nan(length(pdata(p).trials),stimlock_samples);
    pdata(p).y_stimlock=nan(length(pdata(p).trials),stimlock_samples);
        
    %perform preprocessing on each trials
    for tr=1:length(pdata(p).trials) %bp: for each trial in participant p
        %correct sample timing errors
        [x,y,t]=correctSampleTiming(pdata(p).trials(tr).x,pdata(p).trials(tr).y,pdata(p).trials(tr).t);% bp: Correct sampling timings variance to constant sample intervals by linear interpolation
        %mirror right responses to left
        if(pdata(p).response(tr)==2)
            x=screen_width-x;
        end
        %align x coordinates to zero as start point
        x=x-x(1);%bp: this is specified such that the 1st time point x value in each trial is set as 0 (which should be around 640) and everything is relative to that.-
        
        
        %add data to raw data matrices
        pdata(p).x_stimlock(tr,:)=normLength(x,stimlock_samples,0);
        pdata(p).y_stimlock(tr,:)=normLength(y,stimlock_samples,0);
        pdata(p).t_stimlock(tr,:)=normLength(t,stimlock_samples,0);
        %bp: normLength (signal, thelength, interpolate): interpolate or fillup/cut a vector to specific length-
        %-signal=vector or matrix(trial x time) 
        %-thelength: length of data (default:100)
        %-norming procedure:
        %0 => fillup with NaNs; [0,f] => fillup with f; if f==inf, fill up with last value in signal
        %1 => indicates interpolation (default);
    end
     
    %extract all dynamic data from stimulus locked raw data
    [pdata(p).x_warp, pdata(p).y_warp, pdata(p).angle_warp, pdata(p).velocity_warp,...
     pdata(p).dev_warp, pdata(p).angle_stimlock, pdata(p).velocity_stimlock, pdata(p).dev_stimlock]...
        =calcTrajectories(pdata(p).x_stimlock, pdata(p).y_stimlock, warp_samples, sf);
    
    [pdata(p).meandev, pdata(p).maxdev]=calcStatic(pdata(p).x_warp, pdata(p).y_warp);
    %Bp: include the normalized x_stimlock and y_stimlock to calculate meandev: area under the curve and maxdev: maximum deviation from the straight line
end

%% inspect data

% plot data by condition
% aggregate data
clear ca ia pca pia n1r_true n1r_false
% BP: below draws the time course of average movements on X-axis for both conditions

for p=1:length(pdata)
    % OG: average angle warp instead of x_warp --> congruency
    % ca - congruent angle
    % ci - incongruent angle
    ca(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==1,:));%bp: nanmeans extracts mean from a dataset that ignores nan values
    ia(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==2,:));
    
    %OG: average angle_warp --> previous_congruency
    pca(p,:)=nanmean(pdata(p).angle_warp(pdata(p).previous_congruency==1,:));
    pia(p,:)=nanmean(pdata(p).angle_warp(pdata(p).previous_congruency==2,:));
end

% OG: plot mean angles current congruency (OG)
figure;hold on
fig=errorArea(mean(ca),ste(ca),'b');errorArea(mean(ia),ste(ia),'r');
legend('congruent','incongruent');xlabel('warped step');ylabel('angle (radian)');
title('average warp angle for congruent/incongruent (RS)')
ylim([-1, 0.4]);

% OG: plot mean angles previous congruency
fig=figure;hold on
errorArea(mean(pca),ste(pca),'b');errorArea(mean(pia),ste(pia),'r');
legend('prev congruent','prev incongruent');xlabel('warped step');ylabel('angle (radian)');
title('average warp angle for previous congruent/incongruent (RS)')
ylim([-1, 0.4]);

%OG: plot of angle_warp for participant 1 (as an example)

% figure; %plot of warped angles of trials
% plot(pdata(1).angle_warp,1:100);
% title('all warped angle - subject #1')
% xlabel('angle in radian measure')
% ylabel('warp step')

figure; % plot of warped angles for congurent trials
plot(1:100,pdata(1).angle_warp(pdata(1).congruency==1,:)); 
title('congruent warped angle (RS)- subject #1')
ylabel('angle (radian)')
xlabel('warped step')

figure; % plot of warped angles for incongurent trials
plot(1:100,pdata(1).angle_warp(pdata(1).congruency==2,:)); %switched axis: x - angle; y - warp step
title('incongruent warped angle (RS)- subject #1')
ylabel('angle (radian)')
xlabel('warped step')

%OG: sorting and averaging warped angles for each congruency sequence: cc, ic, ci, & ii
clear cc ic ci ii
for p=1:length(pdata)
    n1c=[0;pdata(p).congruency(1:end-1)];
    cc(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==1 & n1c==1,:));%bp: nanmeans extracts mean from a dataset that ignores nan values
    ci(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==2 & n1c==1,:));%bp: nanmeans extracts mean from a dataset that ignores nan values
    ic(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==1 & n1c==2,:));%bp: nanmeans extracts mean from a dataset that ignores nan values
    ii(p,:)=nanmean(pdata(p).angle_warp(pdata(p).congruency==2 & n1c==2,:));%bp: nanmeans extracts mean from a dataset that ignores nan values
end

%OG: plot all four mean angle-arrays in one figure
figure;hold on
errorArea(mean(cc),ste(cc),'b');errorArea(mean(ci),ste(ci),'r');
errorArea(mean(ic),ste(ic),'c');errorArea(mean(ii),ste(ii),'m');
legend('cc','ci','ic','ii');xlabel('warped step');ylabel('angle (radian)');
title('warped angles for cc, ci, ic, & ii (RS)')
ylim([-1, 0.4]);

% check movement quality: how straight is each movement and how many movements show returns (y direction should be increasing constantly (=1)
% BP: Check movement consistency by calculating the movement index (how consistent/straight was the upwards movement, how many backwards movements occurred)
for p=1:length(pdata)
    [cont(p,:),returns(p,:)]=calcMovementContinuity(pdata(p).y_warp,1);
end
figure
boxplot([cont,returns])
set(gca,'XTick',1:2);set(gca,'XTickLabels','continuity|returns')
title('Consistency of movements')

% plot heatmap of pooled trial per condition
% pool all trials
allx=vertcat(pdata.x_warp);
ally=vertcat(pdata.y_warp);
allvel=vertcat(pdata.velocity_warp);
allcong=vertcat(pdata.congruency);

% plot heatmaps of conditions
figure;colormap('hot')

%OG: Changed the scaling --> from e.g. -7:0:2 to -7:0.1:7
% now all deviation trajectory is visible (exceeding the x=2 limitation)
% Also, now the data bins are of size 0.1x0.1 --> the previous step-size of
% "0" didn't work, so matlab presumably put in a default step size, which
% was bigger (-7:0.1:7 now creates an array from -7 to 7 in 0.1 steps)
subplot(1,2,1);
imagep2d(allx(allcong==1,:),-ally(allcong==1,:),-7:0.1:7,-2:0.1:12,[],[],1);
title('congruent')

subplot(1,2,2);
imagep2d(allx(allcong==2,:),-ally(allcong==2,:),-7:0.1:7,-2:0.1:12,[],[],1);
title('incongruent')

% plot heatmap of velocity (=consistency of movement)
figure;colormap('hot')
imagep(allvel,[],[],[],1)
title('velocity (px/ms)')

%% perform TCMR
clear betas
figure
for p = 1:length(pdata)  
    % NOTE BP: Modification are required here depending on number of predictor utilized
    % calc regressors of each participant    
    %OG: congruency effect --> same congruency in trial n and trial n-1
    congs = pdata(p).congruency;% congruency
%     pre_cong = [1;congs(1:end-1)]; % double check, what
%     .previous_congruency is coding
    pre_cong=pdata(p).previous_congruency; % OG: previous congruency
    
    % OG: definition of n1response_x_adaptation regressor
    congs_norm = normScore(congs,[],-1,1);
    pre_cong_norm = normScore(pre_cong,[],-1,1);

    % BP: include 3 way interaction between response bias, previous
    % congruency and current congruency
    cong_x_precong = congs_norm .* pre_cong_norm;
    
    % BP: this creates a data array of 0 and 1 with 1 representing a change in responseN-1 vs responseN, this predictor is used to test if response repetition across trials influence angle trajectory
    % NOTE BP: this should also change depending on predictor numbers
    % concatenate regressors and normalize each to [-1,1]
    
    % OG: n1response excluded as regressor & n1response_x_adaptation included
    regressors = normalizeRegressors([congs, pre_cong, cong_x_precong]); 
    % BP: the code above generates three arrays of regresssor where changes  in responseN-1 vs responseN is coded as 1 and repetition is coded as -1.-
    % -Also congruent trials are coded as -1 and incongruent = 1, so that it
    % matches shift toward target will be positive when we run the analysis.
    % pre_cong is where previous congruent = -1 and previous incongruent = 1
    
    %define data for TCMR
    % OG: Only take the first 80 data points per subject
    % --> reason: presumably "aiming behaviour" at the end of the
    % trajectory (s. figure of warped angles of subject #1)
    % regdata=pdata(p).angle_warp(:,1:80); %limite trajectory data to 80 steps
    regdata=pdata(p).angle_warp;

    %catch plotted information from TCMR
    subplots(length(pdata),p); 
    
    betas(:,:,p)=TCMRegression(regdata,regressors,10); %bp: the 10 is the smoothing variable
    %bp:function[betas,intercept,vifs,fitbetas,fitparams,fitvalue]=TCMR(data,regressors,smoothing,gaussfit,parameterbounds,plotlevel)-
    %-Performs time continuous multiple regression, usually applied to mouse
    %-trajectories, to identify the time profiles of influences (columns of regressors) that were varied
    %-across trials (lines of data and lines of regressors). Smoothing allows
    %-to introduce temporal correlation between time points (comparable to FMRI GLM analyses).
    %-I am guess this code generates a betas matrix that include N (participants) number
    %-of data arraries with N predictors number of row so here the 1st row
    %-is the beta values for response repetition and the 2nd row is congruency
end

% OG: plotting of betas for congruency per subject
figure;
plot(squeeze(betas(1,:,:)));
title('betas for congruency per subject (RS)')
ylabel('beta weight')
xlabel('warped step')


%% determine some process parameters diretly from data
peaks=findStatPeaks(betas,'jackknife');
% bp-detects peaks in the results of TCMRegression
% -Inputs:
%   -betas: the betas from TCMRegression (regressor x time x subject)
%   -detectionmethod:'individual' detects peak of each beta for each subject in betas
%                   'jackknife' detects peak of each beta via jackknifing across subjects
%                   'mean' detects peak of each beta on the average of all subjects 
%
%-Outputs:
%   -peaks: array of structs from peakdetection. The struct contains for each line of betas the fields:
%       -value: the peak height
%      -time: the peak time
%       -sig: the p-value of significance testing against 0
%       -valueSTE: standarderror of value (correct for jackKnifing if necessary)
%       -timeSTE: standarderror of time (correct for jackKnifing if necessary)
%  -peakvalues, peaktimes, peakvaluessig: the raw values for each line of betas for each subject
%               -if 'individual' or for each jackknifed subject if
%               -'jackknife'(use jackKnifeSte, jackKnifeStats, and jackKnifeStats2 
%               -to plot standard errors and perform ttests on these raw data)
segments=findStatSegments(betas,0.05);
% bp:Detects statistically significant segments in the results of TCMRegression
% -Inputs:
%   -betas: the betas from TCMRegression (regressor x time x subject)
%   -alpha: alpha level for significance testing of each time point (ttest against 0)
%
%-Outputs:
%   -segments: cell array of segment start and end times, with one cell for each beta
%   -segmentsraw: matrix (beta  x time) of true/false result of ttest of each time point
%   -segmentstats: matrix (beta  x time) of tvalues from ttest of each time point

%% perform fitting
figure
[fitbetas,fitparams,fitvalues]=fitRegression(betas,'estimate',1);
%bp:Fits gauss curves to betas extracted by TCMR to extract parameters peak time,
%-peak width and peak strength. Peaks in betas are expected to be positive.
%-inputs:
%   -betas = beta weights extracted from TCMRegression, matrix beta x time x subject
%   -paramaterbounds: either matrix with lower & upper parameter boulds gaussian fitting
%                   -parameterbounds:[lower_peaktime,lower_peakwidth,lower_peakstrength;
%                                    upper_peaktime,upper_peakwidth,upper_peakstrength]
%                   -or: String 'default' for default parameter bounds: peaktime [1,signallength], peakwidth [0.05*signallength,signallength], peakstrength [min(betas(:))/1.5,max(betas(:))*1.5]  
%                              'estimate' to estimate parameter bounds from data
%   -plotlevel: 0 = no output
%              -1 = plot aggregated fitted curves for all subjects
%              -2 = plot fitted curves of each subject
%-outputs:
%   -fitbetas: the model betas from gaussian fitting, matrix regressor x time
%   -fitparams: the parameters from gaussian fitting for each regressor: beta x parameter (peak time, peak width (std), peak strength)
%   -fitvalue: R square calculated by mean squared error from gausian fitting procedure normalized to overall variance for each beta

%% save parameters for statistical analyis in external software (e.g. JASP)
saveCSV('fit_parameters.csv',{'n1response_time','n1response_width',...
                              'n1response_strength','cong_time','cong_width','cong_strength'},...
        [squeeze(fitparams(1,:,:))',squeeze(fitparams(2,:,:))']);

%% write table with data to commandline
writeParameterTable(fitparams,fitvalues); %bp: the values in brackets in the output is SEs.
writePeakTable(peaks);
writeSegmentTable(segments);

%% plot results summary for TCMR and fitting results
figure;clear s
s(1)=subplot(1,2,1);
plotRegression(betas);
plotSegmentLines(segments,10);
%bp: segments are identified as significant segements of timepoints identified
plotPeaks(peaks,0.05);
xlabel('time slice');ylabel('\beta weight')
title('RS:Real betas')
s(2)=subplot(1,2,2);
plotRegression(fitbetas);
plotModelLines(fitparams,0.1);
%bp:plotModelLines adds marker lines above a graph (e.g. by plotRegression) indicating the gaussian paramaters extracted by fitRegression
xlabel('time slice');ylabel('\beta weight')
title('RS:Fitted betas')
legend({'Congruency','Previous Congruency', 'Cong x precong'})
linkaxes(s)

