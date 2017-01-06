% Example Script

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

logRet=importdata('Example_Fex.txt');  % load some Data.

dep=logRet(:,1);                    % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep=[constVec logRet(:,2:3)];     % Defining some explanatory variables
k=2;                                % Number of States
S=[0 1 1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.optimizer='fmincon';         % Set opt to fmincon (constrained estimation only works with fmincon)

advOpt.constCoeff.nS_Param{1}={'e'};	

% fix the switching parameter of indep column 2 to 0 at state 1, indep column 3 equal to 0 at state 2 and estimate the rest

advOpt.constCoeff.S_Param{1}={ 0 ,'e' ;	...
                              'e', 0 };		
                          
% estimate both standard deviations/variances

advOpt.constCoeff.covMat{1}(1,1)={0.01^2}; % cells iterate over states
advOpt.constCoeff.covMat{2}(1,1)={'e'};

% estimate only transition probabilities of second state
	
 advOpt.constCoeff.p={0.95,'e'  ; ...
                      0.05,'e' };
                  
[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');
rmpath('data_Files'); 