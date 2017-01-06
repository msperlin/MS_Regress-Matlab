% Likelihood Function for MS(k) Regression

function [sumlik,Output,logLikVec]=MS_Regress_Lik(dep,indep_nS,indep_S,param,k,S,advOpt,disp_out)

useMex=advOpt.useMex;
distrib=advOpt.distrib;
Coeff_Tag=advOpt.Coeff_Tag;
constCoeff=advOpt.constCoeff;

%   Calculation of some preliminary variables

nr=length(dep);
nEq=size(dep,2);

for iEq=1:nEq
    switch distrib
        case 'Normal'
            n_dist_param=1; % Number of distributional parameters
            S_Var{iEq}=S{iEq}(end);   % if model switches in variance
        case 't'
            S_Var{iEq}=S{iEq}(end-1);
            S_df{iEq}=S{iEq}(end);    % if model is switching in degrees of freedom
            n_dist_param=2; % Number of distributional parameters
        case 'GED'
            S_Var{iEq}=S{iEq}(end-1);
            S_k{iEq}=S{iEq}(end);     % if model is switching in k parameter (k as the GED parameter)
            n_dist_param=2; % Number of distributional parameters
    end
    
    n_indep{iEq}=size(indep_nS{iEq},2)+size(indep_S{iEq},2);
    n_S{iEq}=sum(S{iEq}(1:end-n_dist_param));
    n_nS{iEq}=n_indep{iEq}-n_S{iEq};
end

typeCall='estimation';
[Coeff]=param2spec(param,Coeff_Tag,constCoeff,typeCall);

% Making unrestricted parameters positive

switch advOpt.optimizer
    case {'fminsearch','fminunc'}
        for ik=1:k;
            for iEq=1:size(Coeff.covMat{ik},2)
                Coeff.covMat{ik}(iEq,iEq)=abs(Coeff.covMat{ik}(iEq,iEq));
            end
        end
end

% build back the distribution parameters

switch distrib
    case 'Normal'
        myIdx=0;
    case 't'
        myIdx=1;
    case 'GED'
        myIdx=1;
end

for ik=2:k
    for iEq=1:nEq
        if S{iEq}(end-myIdx)==0
            Coeff.covMat{ik}(iEq,iEq)=Coeff.covMat{1}(iEq,iEq); % for cases where variance doesnt switch states
        end
    end
end

% building simetric cov matrix

for ik=1:k
    for iEq=1:nEq
        for jEq=1:nEq
            if iEq<=jEq
                Coeff.covMat{ik}(iEq,jEq)=Coeff.covMat{ik}(jEq,iEq);
            end
        end
    end
end

switch distrib
    case 't'
        for iEq=1:nEq
            for ik=2:k
                if S{iEq}(end)==0
                    Coeff.df{iEq}(1,ik)=Coeff.df{iEq}(1,1); % for cases where variance doesnt switch states
                end
            end
        end
        
    case 'GED'
        for ik=2:k
            for iEq=1:nEq
                if S{iEq}(end)==0
                    Coeff.K{iEq}(1,ik)=Coeff.K{iEq}(1,1); % for cases where variance doesnt switch states
                end
            end
        end
end

for iEq=1:nEq
    
    if n_nS{iEq}==0
        indep_nS{iEq}=zeros(nr,1);
        Coeff.nS_Param{iEq}=0;
    end
    
    if n_S{iEq}==0
        indep_S{iEq}=zeros(nr,1);
        Coeff.S_Param{iEq}=zeros(1,k);
    end
end

% Organizing Coeffs for each state

Cond_mean=cell(nEq,1);
e=cell(nEq,1);
n=zeros(nr,k);

% Vectorized main engine

for i=1:k
    for iEq=1:nEq
        Cond_mean{i}(:,iEq)=indep_nS{iEq}*Coeff.nS_Param{iEq}+indep_S{iEq}*Coeff.S_Param{iEq}(:,i); % Conditional Mean for each equation (cell wise) and each state (column wise)
        e{i}(:,iEq)=dep(:,iEq)-Cond_mean{i}(:,iEq); % Error for each state (cell) for each Equation
    end
    
    switch distrib
        case 'Normal'
            n(:,i)=myMVNPDF(dep,Cond_mean{i},Coeff.covMat{i});
        case 't'

            n(:,i)=( gamma(.5.*(Coeff.df{1}(i)+1)) ) ./ ( (gamma(.5.*Coeff.df{1}(i))).*sqrt(pi().*Coeff.df{1}(i).*Coeff.covMat{i})).* ...
                ((1+((e{i}(:,1)).^2)./(Coeff.df{1}(i).*Coeff.covMat{i})).^(-.5.*(Coeff.df{1}(i)+1)) );  % t density

        case 'GED'
            
            n(:,i)=exp(-1/2.*abs(e{i}(:,1)./sqrt(Coeff.covMat{i})).^(1/Coeff.K{1}(1,i)))./(2.^(Coeff.K{1}(1,i)+1).* ... 
                     sqrt(Coeff.covMat{i}).*gamma(Coeff.K{1}(1,i)+1) );
    end
end

switch advOpt.optimizer
    case {'fminsearch','fminunc'}
        p=ones(k,k);
        temp_p=ones(k,k);
        for i=1:k-1
            for j=1:k
                temp_p(i,j)=normcdf(Coeff.p(i,j));
                p(i+1,j)=p(i,j)*(1-temp_p(i,j));
            end
        end
        p=p.*temp_p;
        Coeff.p=p;
end

if useMex==1
    
    if exist('mex_MS_Filter')==0
        error(['The likelihood function is not being able to use the mex version of the '...
            ' filter. You need to compile the file mex_MS_Filter.cpp in your pc in order for it to work.' ...
            ' More details at pdf document from the zip file.']);
    else
        
        [f,E]=mex_MS_Filter(Coeff.p,[zeros(1,k);n]); 
        f=f(:,2:end)';
        E=E(:,2:end)';
        
    end
else
    
    %     Prealocation of large matrices
    
    E=zeros(nr,k);
    f=zeros(nr,1);
    
    % Setting up first probs of E
    
    firstProb=repmat(1/k,1,k);
    
    for i=1:nr
        
        if i==1 % first probabilities use a naive guess (1/k)
            f(i,1)=ones(k,1)'*(Coeff.p*firstProb'.*n(i,:)'); % MS Filter equation
            E(i,:)=((Coeff.p*firstProb'.*n(i,:)')/f(i,1));
        else
            f(i,1)=ones(k,1)'*(Coeff.p*E(i-1,:)'.*n(i,:)'); % MS Filter equation
            E(i,:)=((Coeff.p*E(i-1,:)'.*n(i,:)')/f(i,1));   % MS Filter equation for probabilities
        end
        
    end
end

% Negative sum of log likelihood for fmincon (fmincon minimizes the function)

sumlik=-sum(log(f(2:end)));

% Control for nan, Inf, imaginary

if isnan(sumlik)||isreal(sumlik)==0||isinf(sumlik)

    % old solution (did not work well)
%     idx1=find(isinf(log(f(2:end)))==1);
%     idx2=find(isnan(log(f(2:end)))==1);
%     
%     idx=[idx1 ;idx2];
%     
%     a=log(f(2:end));
%     a(idx)=[];
%     
%     sumlik=-sum(a);
%     if isempty(a)
%          sumlik=Inf;
%     end

    sumlik=Inf;
    
end

% Building Output structure

Prob_t_1=zeros(nr,k);
Prob_t_1(1,1:k)=1/k; % This is the matrix with probability of s(t)=j conditional on the information in t-1

for i=2:nr
    Prob_t_1(i,1:k)=(Coeff.p*E(i-1,1:k)')'; % prob conditional in t-1
end

f(f==0)=1;

logLikVec=log(f);
Output.Coeff=Coeff;
Output.filtProb=E;
Output.LL=-sumlik;
Output.k=k;
Output.param=param;
Output.S=S;
Output.advOpt=advOpt;

for iEq=1:nEq
    for ik=1:k
        myStdVec{iEq}(:,ik)=repmat(sqrt(Coeff.covMat{ik}(iEq,iEq)),nr,1);
        myCondMean{iEq}(:,ik)=Cond_mean{ik}(:,iEq);
    end
    Output.condMean(:,iEq)=sum(myCondMean{iEq}.*Prob_t_1,2); % conditional mean build with probabiblites conditional in t-1
    Output.condStd(:,iEq)=sum(myStdVec{iEq}.*Prob_t_1,2);
    
    Output.resid(:,iEq)=dep(:,iEq)-Output.condMean(:,iEq);
end

if disp_out==1
    fprintf(1,['\nSum log likelihood for MS Model -> ', num2str(-sumlik)]);
end
