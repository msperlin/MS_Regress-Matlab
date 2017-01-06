nr=size(dep,1);
multIdx=.75;    % controls the factor for starting values of cov matrix
myFactor=1.5;   % controls how to factor each cov matrix for states

% Building bounds for parameters

for ik=1:k
    Coeff0.covMat{ik}=cov(dep).*multIdx;    % initial parameters for cov matrix
    multIdx=multIdx*myFactor;   % for each state, factor cov matrix up
    
    % building bounds for cov matrix
    
    CoeffUpperBnd.covMat{ik}=repmat(Inf,nEq,nEq); 
    CoeffLowerBnd.covMat{ik}=repmat(-Inf,nEq,nEq);
    
    for iEq=1:nEq
        for jEq=1:nEq
            %         CoeffLowerBnd.covMat{i}(j,j)=min(nonzeros(dep(:,j).^2));
            minValue=min( (dep(:,iEq)-mean(dep(:,iEq))).*(dep(:,jEq)-mean(dep(:,jEq))));
            CoeffLowerBnd.covMat{ik}(iEq,jEq)=minValue; % minimum possible covariance value of i and j
            
            maxValue=max( (dep(:,iEq)-mean(dep(:,iEq))).*(dep(:,jEq)-mean(dep(:,jEq))));
            CoeffUpperBnd.covMat{ik}(iEq,jEq)=maxValue; % maximum possible covariance value of i and j
        end
    end
    
end

% building size variables

n_indep=cell(nEq);
n_S=cell(nEq);
n_nS=cell(nEq);
count=cell(nEq);
countS=cell(nEq);

S_S=cell(nEq);
indep_S=cell(nEq);
S_nS=cell(nEq);
indep_nS=cell(nEq);

param0_indep_nS=cell(nEq);
param_ols_S=cell(nEq);
param0_indep_S=cell(nEq);

for iEq=1:nEq
    n_indep{iEq}=size(indep{iEq},2);
    n_S{iEq}=sum(S{iEq}(1:end-n_dist_param));
    n_nS{iEq}=n_indep{iEq}-n_S{iEq};
    count{iEq}=0;
    countS{iEq}=0;
    
    % Checking which parameters will have switching effect
    
    S_S{iEq}=zeros(1,n_S{iEq});
    indep_S{iEq}=zeros(nr,n_S{iEq});
    S_nS{iEq}=zeros(1,n_nS{iEq});
    indep_nS{iEq}=zeros(nr,n_nS{iEq});
    
    for i=1:(length(S{iEq})-n_dist_param)
        
        if S{iEq}(i)==1
            countS{iEq}=countS{iEq}+1;
            S_S{iEq}(countS{iEq})=i;
            indep_S{iEq}(:,countS{iEq})=indep{iEq}(:,i);
        else
            count{iEq}=count{iEq}+1;
            S_nS{iEq}(count{iEq})=i;
            indep_nS{iEq}(:,count{iEq})=indep{iEq}(:,i);
        end
    end
    
    % Calculating starting coefficients (OLS based)
    
    if n_nS{iEq}~=0
        param0_indep_nS{iEq}=regress(dep(:,iEq),indep_nS{iEq}); % simple Ols for param0 of non switching variables
    else
        param0_indep_nS{iEq}=0;
        indep_nS{iEq}=zeros(nr,1);
    end
    
    if n_S{iEq}~=0
        param_ols_S{iEq}=regress(dep(:,iEq),indep_S{iEq}); % simple OlS for param0 of switching variables
        
        param0_indep_S{iEq}=[];
        idx=1;
        for i=0:k-1
            param0_indep_S{iEq}(:,i+1)=idx*param_ols_S{iEq}'; % building param0 of switching variables (changing sign of coefficients)
            idx=idx*-1;
        end
    else
        param0_indep_S{iEq}=0;
        indep_S{iEq}=zeros(nr,1);
    end
    
    % Building the whole param0 as a structure, which will be then translated
    % to vector notation
    
    Coeff0.nS_Param{iEq}=param0_indep_nS{iEq};
    Coeff0.S_Param{iEq}=param0_indep_S{iEq};
    
    switch advOpt.optimizer
        case 'fmincon'
            Coeff0.p=repmat(.1,k,k)+eye(k)*(1-k*.1);
        case {'fminsearch','fminunc'}
            
%             Coeff0.p=repmat(.01,k,k)+eye(k)*(1-k*.01); (OLD VERSION, keep it
%             for future reference)
            
            Coeff0.p=zeros(k-1,k);
            for i=1:k-1
                for j=1:k
                    
                    if i==j
                        Coeff0.p(i,j)=1;
                    else
                        Coeff0.p(i,j)=-1;
                    end
                    
                end
            end
            
    end
    
    for j=1:numel(Coeff0.nS_Param{iEq})   % building max possible values for betas

% OLD CODE for calculation of bounds of parameters (NO LONGER USED)
        
%         mySeries=(dep(:,iEq)-mean(dep(:,iEq))).*(indep_nS{iEq}(:,j)-mean(indep_nS{iEq}(:,j)));
%         myVar=(var(indep_nS{iEq}(:,j)));
%         if myVar==0
%             myVar=1;
%             CoeffUpperBnd.nS_Param{iEq}(j,1)=max(mySeries)./myVar;
%             CoeffLowerBnd.nS_Param{iEq}(j,1)=min(mySeries)./myVar;
%         else
%             CoeffUpperBnd.nS_Param{iEq}(j,1)=max(mySeries)./myVar;
%             CoeffLowerBnd.nS_Param{iEq}(j,1)=min(mySeries)./myVar;
%         end

% Parameters are free

        CoeffUpperBnd.nS_Param{iEq}(j,1)=inf;
        CoeffLowerBnd.nS_Param{iEq}(j,1)=-inf;
    end
    
    for j1=1:size(Coeff0.S_Param{iEq},1)   % building max possible values for betas
        for j2=1:size(Coeff0.S_Param{iEq},2)
            
% OLD CODE for calculation of parameter's boundaries (NO LONGER USED)
            
%             mySeries=(dep(:,iEq)-mean(dep(:,iEq))).*(indep_S{iEq}(:,j1)-mean(indep_S{iEq}(:,j1)));
%             myVar=(var(indep_S{iEq}(:,j1)));
%             if myVar==0 % case where indep variable is a constant
%                 myVar=1;
%                 CoeffUpperBnd.S_Param{iEq}(j1,j2)=max(dep(:,iEq));
%                 CoeffLowerBnd.S_Param{iEq}(j1,j2)=min(dep(:,iEq));
%             else
%                 CoeffUpperBnd.S_Param{iEq}(j1,j2)=max(mySeries)./myVar;
%                 CoeffLowerBnd.S_Param{iEq}(j1,j2)=min(mySeries)./myVar;
%             end

% parameters are fre
                CoeffUpperBnd.S_Param{iEq}(j1,j2)=inf;
                CoeffLowerBnd.S_Param{iEq}(j1,j2)=-inf;
        end
    end
    
    switch advOpt.optimizer
        case 'fmincon'
            CoeffUpperBnd.p=ones(k,k); % Fixed (value k-1 for new algorithm)
            CoeffLowerBnd.p=zeros(k,k);
        case {'fminsearch','fminunc'}
            CoeffUpperBnd.p=ones(k-1,k); % Fixed (value k-1 for new algorithm)
            CoeffLowerBnd.p=zeros(k-1,k);
    end
    
    switch distrib
        
        case 't'
            Coeff0.df{iEq}=repmat(10,1,k);
            CoeffUpperBnd.df{iEq}=Inf(1,k);
            CoeffLowerBnd.df{iEq}=repmat(0.1,1,k);
        case 'GED'
            Coeff0.K{iEq}=repmat(.5,1,k);
            CoeffUpperBnd.K{iEq}=repmat(5,1,k);
            CoeffLowerBnd.K{iEq}=repmat(0.01,1,k);
    end
    
    
end

[Coeff_Tag,param0]=spec2param(Coeff0);      % converting starting coefficients to vector notation
[~,uB]=spec2param(CoeffUpperBnd);   % converting upper bound to vector notation
[~,lB]=spec2param(CoeffLowerBnd);   % converting lower bound to vector notation

%     procedures for adjusting parameter vector (for estimated/non estimated
%     coefficients)

% building a new tag structure for the coefficients

newCoeff_Tag.p=zeros(size(Coeff_Tag.p));
for iEq=1:nEq
    
    
    newCoeff_Tag.nS_Param{iEq}=zeros(size(Coeff_Tag.nS_Param{iEq}));
    newCoeff_Tag.S_Param{iEq}=zeros(size(Coeff_Tag.S_Param{iEq}));
    
    switch distrib
        case 't'
            newCoeff_Tag.df{iEq}=zeros(size(Coeff_Tag.df{iEq}));
        case 'GED'
            newCoeff_Tag.K{iEq}=zeros(size(Coeff_Tag.K{iEq}));
    end
end

for ik=1:k
    newCoeff_Tag.covMat{ik}=zeros(size(Coeff_Tag.covMat{ik}));
end

% calculating new tags

idxVec=zeros(0);
count_e=1;
count_ne=1;

for ik=1:k
    for i=1:nEq
        for j=1:nEq
            if ~(isnumeric(constCoeff.covMat{ik}{i,j}))
                newCoeff_Tag.covMat{ik}(i,j)=count_e;
                count_e=count_e+1;
            else
                idxVec(count_ne)=Coeff_Tag.covMat{ik}(i,j);
                count_ne=count_ne+1;
            end
        end
    end
end

for iEq=1:nEq
    for i=1:numel(Coeff0.nS_Param{iEq})
        if ~(isnumeric(constCoeff.nS_Param{iEq}{i}))
            newCoeff_Tag.nS_Param{iEq}(i,1)=count_e;
            count_e=count_e+1;
        else
            idxVec(count_ne)=Coeff_Tag.nS_Param{iEq}(i);
            count_ne=count_ne+1;
        end
    end
end

for iEq=1:nEq
    for i=1:size(Coeff0.S_Param{iEq},1)
        for j=1:size(Coeff0.S_Param{iEq},2)
            if ~(isnumeric(constCoeff.S_Param{iEq}{i,j}))
                newCoeff_Tag.S_Param{iEq}(i,j)=count_e;
                count_e=count_e+1;
            else
                idxVec(count_ne)=Coeff_Tag.S_Param{iEq}(i,j);
                count_ne=+count_ne+1;
            end
        end
    end
end

for i=1:size(Coeff0.p,1)
    for j=1:size(Coeff0.p,2)
        if ~(isnumeric(constCoeff.p{i,j}))
            newCoeff_Tag.p(i,j)=count_e;
            count_e=count_e+1;
        else
            idxVec(count_ne)=Coeff_Tag.p(i,j);
            count_ne=count_ne+1;
        end
    end
end

for iEq=1:nEq
    switch distrib
        case 't'
            for i=1:size(Coeff0.df{iEq},2)
                if ~(isnumeric(constCoeff.df{iEq}{i}))
                    newCoeff_Tag.df{iEq}(i)=count_e;
                    count_e=count_e+1;
                else
                    idxVec(count_ne)=Coeff_Tag.df{iEq}(i);
                    count_ne=count_ne+1;
                end
            end
            
        case 'GED'
            
            for i=1:size(Coeff0.K{iEq},2)
                if ~(isnumeric(constCoeff.K{iEq}{i}))
                    newCoeff_Tag.K{iEq}(i)=count_e;
                    count_e=count_e+1;
                else
                    idxVec(count_ne)=Coeff_Tag.K{iEq}(i);
                    count_ne=count_ne+1;
                end
            end
    end
end

% Cleaning values in vectors that won't be estimated

param0(idxVec)=[];
lB(idxVec)=[];
uB(idxVec)=[];

Coeff_Tag=newCoeff_Tag;
advOpt.Coeff_Tag=Coeff_Tag;