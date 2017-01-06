% Function for forecasting in t+1 an Markov Switching regression model
% estimated with MS_Regress_Fit.m
%
%   Input:  Spec_Output - Specification output from estimation (check MS_Regress_Fit.m)
%           newIndepData - New information that has arrived for t+1 (maybe lagged variables ?)
%
%   Output: meanFor - Forecast for the mean equation (column iterating over equations of system)
%           stdFor - Forecast for the standard deviation (columns iterating over
%           equations of system)
% 
%   Author: Marcelo Perlin
%   Email:  marceloperlin@gmail.com

function [meanFor,stdFor]=MS_Regress_For(Spec_Out,newIndepData)

% Retrieving variables from Spec_Output

nEq=size(Spec_Out.condMean,2);
k=Spec_Out.k;
S=Spec_Out.S;
Coeff=Spec_Out.Coeff;

for ik=1:k
    myVariance{ik}=(diag(Coeff.covMat{ik}));
end

distrib=Spec_Out.advOpt.distrib;

for iEq=1:nEq
    switch distrib
        case 'Normal'
            S_Std{iEq}=S{iEq}(end);
            n_dist_param=1; % Number of d
        case 't'
            S_Std{iEq}=S{iEq}(end-1);
            n_dist_param=2; % Number of distributional parameters
        case 'GED'
            S_Std{iEq}=S{iEq}(end-1);
            n_dist_param=2; % Number of distributional parameters
    end
end
filtProb=Spec_Out.filtProb;
p=Spec_Out.Coeff.p;

indep_S=cell(nEq,1);
indep_nS=cell(nEq,1);
for iEq=1:nEq
    
    count_S=0;
    count_nS=0;
    for i=1:length(S{iEq})-n_dist_param
        
        if S{iEq}(i)==1
            count_S=count_S+1;
            indep_S{iEq}(:,count_S)=newIndepData(:,i);
        else
            count_nS=count_nS+1;
            indep_nS{iEq}(:,count_nS)=newIndepData(:,i);
        end
    end
    
    if count_nS==0
        indep_nS{iEq}=0;
    end
    
end


% Building Forecasts

for iEq=1:nEq
    for j=1:k
        
        meanFor_S{iEq}(1,j)=indep_nS{iEq}*Coeff.nS_Param{iEq}+indep_S{iEq}*Coeff.S_Param{iEq}(:,j); % mean forecast for each state
        myVariance_S{iEq}(1,j)=myVariance{j}(iEq); % sigma^2 forecast for each state
        
    end
    
    meanFor(1,iEq)=meanFor_S{iEq}*(p*filtProb(end,:)'); % mean t+1 forecast
    stdFor(1,iEq)=sqrt(myVariance_S{iEq}*(p*filtProb(end,:)'));   % std t+1 forecast

end


