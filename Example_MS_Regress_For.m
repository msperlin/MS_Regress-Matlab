% Example Script using the function MS_Regress_For

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

logRet=importdata('Example_Fex.txt');  % load some Data.

idx=501;    % this is an index of where to forecast the new observation ( in this case the obs 500)

dep=logRet(1:idx-1,1);              % Load only data up to idx-1
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep=[constVec logRet(1:idx-1,2:3)];% Defining some explanatory variables (only info up to idx-1)
k=2;                                % Number of States
S=[1 0 0 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

newIndepData=[1 logRet(idx,2:3)]; % this is the new information feeded to the forecasting function

[meanFor,stdFor]=MS_Regress_For(Spec_Out,newIndepData);

% printing results to screen

fprintf(1,'\nThe mean forecast (t+1 conditional on t) for obs %i is %4.4f',idx,meanFor);
fprintf(1,'\nThe sigma forecast (t+1 conditional on t) for obs %i is %4.4f\n',idx,stdFor);

rmpath('m_Files');
rmpath('data_Files');
