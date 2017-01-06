% Script for simulating and fitting a 2 state MS regression model with two equations. 
% Just press f5 to run it

addpath('m_Files'); % add 'm_Files' folder to the search path

clear;

nr=1000;            % Number of observations in simulation
distrib='Normal';   % The distribution assumption ('Normal' or 't')

Coeff.p=[.9 .1  ; ...    % Transition matrix (this also defines the value of k)
         .1 .9 ];

% Setting up Coefficients for Equation # 1 of system

Coeff.S{1}=[1 0 0];  % Setting up which variables at indep will have switching effect

Coeff.nS_param{1}(1,1)= .2;    % Setting up the coefficients at non switching parameters
Coeff.nS_param{1}(2,1)=-.2;    % Setting up the coefficients at non switching parameters

Coeff.S_param{1}(:,1)= .5;    % Setting up the coefficients at switching parameters
Coeff.S_param{1}(:,2)=-.3;    % Setting up the coefficients at switching parameters

Coeff.indepMean{1}=[0 0 0];
Coeff.indepStd{1}=[.1 .1 .1];

% Setting up Coefficients for Equation # 2 of system

Coeff.S{2}=[1 0 0];  % Setting up which variables at indep will have switching effect

Coeff.nS_param{2}(1,1)= .5;    % Setting up the coefficients at non switching parameters
Coeff.nS_param{2}(2,1)=-.5;    % Setting up the coefficients at non switching parameters

Coeff.S_param{2}(:,1)= .3;    % Setting up the coefficients at switching parameters
Coeff.S_param{2}(:,2)=-.7;    % Setting up the coefficients at switching parameters

Coeff.indepMean{2}=[0 0 0];
Coeff.indepStd{2}=[.1 .1 .1];

% Setting up CoVariance MAtrix of system, where cell iterates over states
% (and not equations)

Coeff.covMat{1}=[0.7697    0.3081 ; ...       % a positive definite matrix
                0.3081    0.4989]*1E-3;

Coeff.covMat{2}=[0.7697    -0.3081 ; ...
                 -0.3081    0.4989]*1E-3;

k=size(Coeff.p,1);  % getting the value of k, according to Coeff.p

[Simul_Out]=MS_Regress_Sim(nr,Coeff,k,distrib);

dep=Simul_Out.dep;               % Defining dependent variable from .mat file
indep{1}=Simul_Out.indep{1};     % Defining some explanatory variables
indep{2}=Simul_Out.indep{2};     % Defining some explanatory variables

k=2;                                % Number of States
S{1}=[1 0 0 1];                     % Defining which parts of the equation #1 will switch states
S{2}=[1 0 0 1];                     % Defining which parts of the equation #2 will switch states

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;                % estimate whole cov matrix

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');