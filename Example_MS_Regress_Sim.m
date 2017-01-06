% Example Script for MS_Regress_Simul.m

addpath('m_Files'); % add 'm_Files' folder to the search path

clear; clc;

nr=500;                    % Number of observations in simulation
advOpt.distrib='Normal';   % The distribution assumption ('Normal' or 't')

% Transition matrix (this also defines the value of k)
Coeff.p=[.95 .1; ...    
         .05 .9];

% Setting up which variables at indep will have switching effect
Coeff.S=[1 0];  

% Setting up the coefficients at non switching parameters (each row is each
% variables coefficient). The order is the same as Coeff.S

% Setting up the coefficients at non switching parameters 
Coeff.nS_param=0;    
 
% Setting up the coefficients at non switching parameters (each row is each
% variables coefficient and each collum is each state). This example has
% 1 switching parameter and 2 states

Coeff.S_param(1,1)= .5;    
Coeff.S_param(1,2)=-.5;    

% Setting up the standard deviavion of the model at each state 
Coeff.Std(1,1)=0.5;  
Coeff.Std(1,2)=1;  

% The explanatory variables used in the simulation are always random normal, with
% specific mean and standard deviation 

Coeff.indepMean=[1 0];  
Coeff.indepStd= [0 0];  

% getting the value of k, according to Coeff.p
k=size(Coeff.p,1);  

[Simul_Out]=MS_Regress_Sim(nr,Coeff,k,advOpt.distrib); % calling simulation function

rmpath('m_Files');