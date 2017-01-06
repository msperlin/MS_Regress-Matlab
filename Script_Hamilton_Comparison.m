% Script for a 2 step estimation of Hamilton (1989) model 
% Data available from:
% http://weber.ucsd.edu/~jhamilto/software.htm

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files'); % add 'm_Files' folder to the search path

gnpVec=importdata('gnp_Hamilton.txt');    % load gnp Data

gnpGrowth=(gnpVec(2:end)-gnpVec(1:end-1))./gnpVec(1:end-1);

dep=gnpGrowth;

constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep=constVec;                     % Defining regressors

k=2;                                % Number of States
S=[1 0];                            % Defining which parts of the equation will switch states (column 1 and variance only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

resid=dep-Spec_Out.condMean;    % Retrieving residuals

[myOLS se r]=regress(resid(5:end),[resid(4:end-1),resid(3:end-2),resid(2:end-3),resid(1:end-4)]);   % using OLS for estimating phi vector

myStd=std(r);   % this standard deviation aproximates the standard deviation of Hamilton's model

rmpath('m_Files');
rmpath('data_Files'); % add 'm_Files' folder to the search path