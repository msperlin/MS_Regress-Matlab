% Example Script for MS_Regress_Fit.m with mex version of Hamilton's
% Filter. 
%
% The mex function was created using a c++ API provided by NR
% (http://www.nr.com/). This makes it easier  by using matrix notation (and
% not pointers notation in the usual C enviroment).
%
% The mex function was developed using matlab 2008a and C++ MS Visual Express 2008 compiler
% (you can get it free at microsoft site). WARNING. It will NOT compile under LCC. 

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

logRet=importdata('Example_Fex.txt');  % load some Data.

dep=logRet(:,1);                    % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep=[constVec logRet(:,2:3)];     % Defining some explanatory variables
k=2;                                % Number of States
S=[1 0 0 1];                        % Defining which ones from indep will have switching effect (in this case variable 1 (constant), only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors
advOpt.useMex=1;                    % Defining use of mex version of likelihood function 

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');
rmpath('data_Files');

