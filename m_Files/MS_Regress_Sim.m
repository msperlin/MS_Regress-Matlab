% Function for simulation of an Markov Switching regression model
%
%   Input:  nr      - Number of rows (number of time periods to simulate)
%           Coeff   - Coeff structure with all the coefficients
%                       More details of how to build the Coeff struct at:
%                       Example_Script_MS_Regress_Simul_2_States.m
%           k       - Number of States
%           distrib - The distribution assumption - 'Normal' or 't'
%
%   Output: Simul_Out - A structure with following fields:
%
%               dep - Simulated explained time series
%               indep - the explanatory series
%               States- The "true" simulated states through time
%               Coeff - A structure with all coefficients
%
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   PhD Student in finance ICMA/UK
%   Created: August/2007

function [Simul_Out]=MS_Regress_Sim(nr,Coeff,k,distrib)

if nargin()==3
    distrib='Normal';
end

if strcmp(distrib,'Normal')==0&&strcmp(distrib,'t')==0
    error('The distrib input should be ''Normal'' or ''t''');
end

if nr<0
    error(['Wow, a negative number of observations! Are you implying that you live in a different dimension where time units decreases instead of increasing, and that we ' ...
        ' are in a convergence period (t=t*). If you tell me all about the future, I''ll tell you all about the past.']);
end

if any( (sum(Coeff.p)>1.0001)|(sum(Coeff.p)<.999) )
    error('The sum of each collum in Coeff.p should be equal to 1 (they represent probabilities of a full process.)')
end

if any([Coeff.p<0 Coeff.p>1])
    error('The Coeff.p stores probabilities and they should be lower than 1 and higher than 0, unless you have a new theory about probabilities (it should be a very convincing one)');
end

if iscell(Coeff.S)
    multiVarFlag=1;
    nEq=numel(Coeff.S);
else
    multiVarFlag=0;
end



if ~multiVarFlag
    
    if ~isfield(Coeff,'nS_param')
        Coeff.nS_param=[];
    end
            
    if size(Coeff.S,2)~=size(Coeff.indepMean,2)
        error('The Coeff.S has a different number of elements than the Coeff.indepMean.');
    end

    if size(Coeff.indepStd,2)~=size(Coeff.indepMean,2)
        error('The Coeff.indepStd has a different number of collum than the value of k.');
    end

    if any(Coeff.Std<0)||any(Coeff.indepStd<0)
        error('All values at Coeff.Std and Coeff.indepStd should be positive (they are standard deviations)');
    end

    if Coeff.S~=0&Coeff.S~=1
        error('The Coeff.S should only have values 1 and 0');
    end

    if size(Coeff.nS_param,1)~=sum(Coeff.S==0)
        error('The Coeff.nS_param should have the same nubmer of rows as the number of zero elements at Coeff.S ');
    end

    if size(Coeff.S_param,1)~=sum(Coeff.S==1)
        error('The Coeff.S_param should have the same number of rows as number of 1 elements at Coeff.S ');
    end

    if size(Coeff.S_param,2)~=k
        error('The Coeff.S_param should have the same number of collums as k ');
    end

    if isempty(Coeff.nS_param)

        Coeff.nS_param=0;

    end

else

    for i=1:nEq
        if ~isfield(Coeff,'nS_param')
            Coeff.nS_param{i}=[];
        end
        
        if size(Coeff.S{i},2)~=size(Coeff.indepMean{i},2)
            error(sprintf('Coeff.S{%i} has a different number of elements than the Coeff.indepMean.',i));
        end

        if size(Coeff.indepStd{i},2)~=size(Coeff.indepMean{i},2)
            error(sprintf('The Coeff.indepStd{%i} has a different number of collum than the value of k.',i));
        end

        if (Coeff.S{i}~=0)&(Coeff.S{i}~=1)
            error(sprintf('The Coeff.S{%i} should only have values 1 and 0',i));
        end

        if size(Coeff.nS_param{i},1)~=sum(Coeff.S{i}==0)
            error(sprintf('The Coeff.nS_param{%i} should have the same nubmer of rows as the number of zero elements at Coeff.S ',i));
        end

        if size(Coeff.S_param{i},1)~=sum(Coeff.S{i}==1)
            error(sprintf('Coeff.S_param{%i} should have the same number of rows as number of 1 elements at Coeff.S ',i));
        end

        if size(Coeff.S_param{i},2)~=k
            error(sprintf('The Coeff.S_param{%i} should have the same number of collums as k ',i));
        end

        if isempty(Coeff.nS_param{i})
            Coeff.nS_param{i}=0;
        end


    end
end

Rnd=rand(nr,1); % random seed for states transition (uniform probabilities)

% Prealocation of large matrices

States=zeros(nr,k);

States(1,1)=1;  % starting with state 1 for first obs

for i=2:nr  % Loop for construction of states (maybe vectorize later ??)

    state_past=find(States(i-1,:)==1);

    if Rnd(i,1)<Coeff.p(state_past,state_past)

        States(i,state_past)=1; % when staying at last state

    else    % when changing to other states

        idx_other=find(States(i-1,:)==0);
        Prob2=Coeff.p(:,state_past);

        a=[Coeff.p(state_past,state_past) ; Prob2(idx_other)];

        cum_sum=cumsum(a);
        sorted=sort([cum_sum ; Rnd(i,1)]); % throw the prob at cumsum of other states to get
        % where it stands (where to
        % switch)

        idx=find(Rnd(i,1)==sorted)-1;      % find index

        States(i,idx_other(idx))=1;        % change state

    end
end

% Creation of random explanatory variables

if multiVarFlag

    for i=1:nEq
        for j=1:length(Coeff.S{i})
            indep{i}(:,j)=Coeff.indepMean{i}(j)+Coeff.indepStd{i}(j)*randn(nr,1);
        end
    end

else

    for i=1:length(Coeff.S)
        indep(:,i)=Coeff.indepMean(i)+Coeff.indepStd(i)*randn(nr,1);
    end
end

% Checking which parameters will have switching effect

if multiVarFlag

    for i=1:nEq
        count{i}=0;
        countS{i}=0;
    end

    for iEq=1:nEq

        for i=1:length(Coeff.S{iEq})
            if Coeff.S{iEq}(i)==1
                countS{iEq}=countS{iEq}+1;
                indep_S{iEq}(:,countS{iEq})=indep{iEq}(:,i);
            else
                count{iEq}=count{iEq}+1;
                indep_nS{iEq}(:,count{iEq})=indep{iEq}(:,i);
            end
        end
    end

else
    count=0;
    countS=0;

    for i=1:length(Coeff.S)

        if Coeff.S(i)==1
            countS=countS+1;
            indep_S(:,countS)=indep(:,i);
        else
            count=count+1;
            indep_nS(:,count)=indep(:,i);
        end
    end
end

% creating model's error for each state

if multiVarFlag
    for i=1:k
        rndError{i}=mvnrnd(zeros(1,nEq),Coeff.covMat{i},nr);
    end
else

    for i=1:k
        switch distrib
            case 'Normal'
                rndError(:,i)=Coeff.Std(i).*randn(nr,1);
            case 't'
                rndError(:,i)=Coeff.Std(i).*trnd(Coeff.v(1,i),nr,1);
        end
    end
end
% conditional mean (the simulated series) calcualation

if multiVarFlag

    for i=1:nEq
        for j=1:k
            tempRndError{i}(:,j)=rndError{j}(:,i);
        end
    end

    for i=1:nEq
        Sim_x(:,i)=indep_nS{i}*Coeff.nS_param{i}+sum(States.*(indep_S{i}*Coeff.S_param{i}),2)+sum(States.*tempRndError{i},2);
    end
else
    Sim_x=indep_nS*Coeff.nS_param+sum(States.*(indep_S*Coeff.S_param),2)+sum(States.*rndError,2);
end

% Passing up to output structure

Simul_Out.dep=Sim_x;
Simul_Out.Coeff=Coeff;
Simul_Out.States=States;
Simul_Out.indep=indep;
Simul_Out.k=k;
Simul_Out.S=Coeff.S;

% Plotting simulated series

if multiVarFlag
    for i=1:nEq
        myLeg{i}=['Simulated MS Series #' num2str(i)];
    end
else
    myLeg='Simulated MS Series';
end

figure(1);
plot(Simul_Out.dep);
legend(myLeg);
xlabel('Time');
ylabel('Simulated MS Series');