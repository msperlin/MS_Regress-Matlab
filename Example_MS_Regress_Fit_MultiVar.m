% Example Script MS_Regress_Fit.m

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

logRet=importdata('Example_Fex.txt');  % load some Data.

dep=logRet(:,1:2);                  % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep{1}=constVec;                  % Defining some explanatory variables
indep{2}=constVec;                  % Defining some explanatory variables
k=2;                                % Number of States
S{1}=[1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
S{2}=[1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');
rmpath('data_Files'); 