% Function for estimation of a general Markov Switching regression
%
%   Input:  dep     - Dependent Variable (vector (univariate model) or matrix (multivariate) )
%           indep   - Independent variables (explanatory variables), should
%                     be cell array in the case of multivariate model (see examples).
%           k       - Number of states (integer higher or equal to 2)
%           S       - This variable controls for where to include a Markov Switching effect.
%                     See pdf file for details.
%           advOpt  - A structure with advanced options for algorithm.
%                     See pdf file for details.
%
%   Output: Spec_Output - A structure with all information regarding the
%                         model estimated from the data (see pdf for details).
%
%   Author: Marcelo Perlin (UFRGS/BR)
%   Contact:  marceloperlin@gmail.com

function [Spec_Output]=MS_Regress_Fit(dep,indep,k,S,advOpt)

% Error checking lines

checkInputs(); % checking if inputs variables are OK

% building constCoeff for the cases when it is not specified

build_constCoeff();

% checking if all fields are specified and make sense

check_constCoeff();

% checking sizes of fields in constCoeff

checkSize_constCoeff();

% Pre calculations before calling the optimizer

preCalc_MSModel();

% Initialization of optimization algorithm

warning('off');

options=optimset('fmincon');
options=optimset(options,'display','off');

dispOut=advOpt.printIter;

% Defining linear contraints in model

A=[];   % inequality constrain (not used)
b=[];   % inequality constrain (not used)

% equality constraint (each collum of Coeff.p must sum to 1)

beq=ones(k,1);  
Aeq=zeros(k,numel(param0));

for i=1:k
    idx=Coeff_Tag.p(:,i);
    
    for j=1:numel(idx)
        if idx(j)==0
            continue;
        else
            Aeq(i,idx(j))=1;
        end
    end
    
end

for i=1:k
    if all(Aeq(i,:)==0)
        Aeq(i,:)=0;    % fixing equality restrictions for when using contrained estimation in Coeff.p
        beq(i,:)=0;
    end
end

param0=param0'; % changing notation for param0

% Call to optimization function

switch advOpt.optimizer
    case 'fminsearch'
        options=optimset('fminsearch');
        options=optimset(options,'display','off');
        options=optimset(options,'MaxIter',500*numel(param0));
        options=optimset(options,'MaxFunEvals',500*numel(param0));
        
        [param]=fminsearch(@(param)MS_Regress_Lik(dep,indep_nS,indep_S,param,k,S,advOpt,dispOut),param0,options);
        
    case 'fminunc'
        options=optimset('fminunc');
        options=optimset(options,'display','off');
        [param]=fminunc(@(param)MS_Regress_Lik(dep,indep_nS,indep_S,param,k,S,advOpt,dispOut),param0,options);
       
    case 'fmincon'
        options=optimset('fmincon');
        options=optimset(options,'display','off');
        [param]=fmincon(@(param)MS_Regress_Lik(dep,indep_nS,indep_S,param,k,S,advOpt,dispOut),param0, ...
            A,b,Aeq,beq,lB,uB,[],options);
        
        
end

% Calculation of Covariance Matrix

[V]=getvarMatrix_MS_Regress(dep,indep_nS,indep_S,param,k,S,std_method,advOpt);
param_std=sqrt(diag((V)));

% Controls for covariance matrix. If found imaginary number for variance, replace with
% Inf. This will then be showed at output

param_std(isinf(param_std))=0;
param_pvalues=2*(1-tcdf(abs(param./param_std),nr-numel(param)));

if ~isreal(param_std)
    for i=1:numel(param)
        if ~isreal(param_std(i))
            param_std(i)=Inf;
        end
    end
end

typeCall='se_calculation';

[Coeff_SE]=param2spec(param_std,Coeff_Tag,constCoeff,typeCall);
[Coeff_pValues]=param2spec(param_pvalues,Coeff_Tag,constCoeff,typeCall);

% After finding param, filter it to the data to get estimated output

[sumlik,Spec_Output]=MS_Regress_Lik(dep,indep_nS,indep_S,param,k,S,advOpt,0);

% calculating smoothed probabilities

Prob_t_1=zeros(nr,k);
Prob_t_1(1,1:k)=1/k; % This is the matrix with probability of s(t)=j conditional on the information in t-1

for i=2:nr
    Prob_t_1(i,1:k)=(Spec_Output.Coeff.p*Spec_Output.filtProb(i-1,1:k)')';
end

filtProb=Spec_Output.filtProb;

P=abs(Spec_Output.Coeff.p);

smoothProb=zeros(nr,k);
smoothProb(nr,1:k)=Spec_Output.filtProb(nr,:);  % last observation for starting filter

for i=nr-1:-1:1     % work backwards in time for smoothed probs
    for j1=1:k
        for j2=1:k
            smooth_value(1,j2)=smoothProb(i+1,j2)*filtProb(i,j1)*P(j2,j1)/Prob_t_1(i+1,j2);
        end
        smoothProb(i,j1)=sum(smooth_value);
    end
end

% Calculating Expected Duration of regimes

stateDur=1./(1-diag(Spec_Output.Coeff.p));
Spec_Output.stateDur=stateDur;

% passing values to output structure

Spec_Output.smoothProb=smoothProb;
Spec_Output.nObs=size(Spec_Output.filtProb,1);
Spec_Output.nEq=nEq;
Spec_Output.Number_Parameters=numel(param);
Spec_Output.advOpt.distrib=distrib;
Spec_Output.advOpt.std_method=std_method;
Spec_Output.Coeff_SE=Coeff_SE;
Spec_Output.Coeff_pValues=Coeff_pValues;
Spec_Output.AIC=2*numel(param)-2*Spec_Output.LL;
Spec_Output.BIC=-2*Spec_Output.LL+numel(param)*log(Spec_Output.nObs*nEq);

% ploting probabilities

if advOpt.doPlots
    doPlots();
end

% Sending output to matlab's screen

disp(' ');
if advOpt.printOut
    doOutScreen()
end
