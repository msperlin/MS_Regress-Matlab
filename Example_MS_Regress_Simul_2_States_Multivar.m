% Script for simulating a 2 state MS regression model with 2 equations. 
% Just press f5 to run it

addpath('m_Files'); % add 'm_Files' folder to the search path

clear; 

nr=1000;            % Number of observations in simulation
distrib='Normal';   % The distribution assumption ('Normal' or 't')

Coeff.p=[.8 .1 ; ...    % Transition matrix (this also defines the value of k)
         .2 .9 ];

% Setting up Coefficients for Equation # 1 of system     
     
Coeff.S{1}=[1 0 0];  % Setting up which variables at indep will have switching effect

Coeff.nS_param{1,1}(1,1)= .2;    % Setting up the coefficients at non switching parameters 
Coeff.nS_param{1,1}(2,1)=-.2;    % Setting up the coefficients at non switching parameters 

Coeff.S_param{1,1}(:,1)= .5;    % Setting up the coefficients at switching parameters 
Coeff.S_param{1,1}(:,2)=-.3;    % Setting up the coefficients at switching parameters 

Coeff.indepMean{1,1}=[0 0 0];   % mean of independent variables in eq 1
Coeff.indepStd{1,1}=[1 1 1];    % std of independent variables in eq 1

% Setting up Coefficients for Equation # 2 of system     
     
Coeff.S{2}=[1 0 0];  % Setting up which variables at indep will have switching effect

Coeff.nS_param{2,1}(1,1)= .5;    % Setting up the coefficients at non switching parameters 
Coeff.nS_param{2,1}(2,1)=-.5;    % Setting up the coefficients at non switching parameters 

Coeff.S_param{2,1}(:,1)= .1;    % Setting up the coefficients at switching parameters 
Coeff.S_param{2,1}(:,2)=-.1;    % Setting up the coefficients at switching parameters 

Coeff.indepMean{2,1}=[0 0 0];   % mean of independent variables in eq 2
Coeff.indepStd{2,1}=[1 1 1];    % std of independent variables in eq 2

% Setting up CoVariance MAtrix of system (iterates over states)

Coeff.covMat{1,1}=[5 -4 ; ...
                   -4 5];

Coeff.covMat{1,2}=[5 2 ; ...
                   2 5];             
             
k=size(Coeff.p,1);  % getting the value of k, according to Coeff.p

[Simul_Out]=MS_Regress_Sim(nr,Coeff,k,distrib);

rmpath('m_Files');