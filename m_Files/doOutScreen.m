
if nEq>1
    typeModel='Multivariate';
else
    typeModel='Univariate';
end

fprintf(1,'\n\n***** Numerical Optimization Converged *****\n\n');
fprintf(1,['Final log Likelihood: ',num2str(Spec_Output.LL),'\n']);
fprintf(1,['Number of estimated parameters: ',num2str(Spec_Output.Number_Parameters),'\n']);
fprintf(1,['Number of Observations: ',num2str(Spec_Output.nObs),'\n']);
fprintf(1,['Number of Equations: ',num2str(Spec_Output.nEq),'\n']);
fprintf(1,['Optimizer: ',Spec_Output.advOpt.optimizer,'\n']);
fprintf(1,['Type of Switching Model: ',typeModel,'\n']);
if mvFlag
    fprintf(1,['Number of Equations in System: ',num2str(nEq),'\n']);
end
fprintf(1,['Distribution Assumption -> ',distrib,'\n']);
fprintf(1,['Method SE calculation -> ',num2str(Spec_Output.advOpt.std_method)]);

for iEq=1:nEq
    fprintf(1,['\n\n***** Final Parameters for Equation #' num2str(iEq) ' ***** \n']);
    
    fprintf(1,'\n---> Non Switching Parameters <---\n');
    
    if (n_nS{iEq}==0)&&(any([S_Var{iEq},S_df{iEq},S_K{iEq}])==0)
        fprintf(1,'\nThere was no Non Switching Parameters for Indep matrix of Equation #%i. Skipping this result',iEq);
    else
        for i=1:n_nS{iEq}
            fprintf(1,'\nNon Switching Parameter for Equation #%i, Indep column %i ',iEq, S_nS{iEq}(i));
            fprintf(1,['\n     Value:                ', num2str(Spec_Output.Coeff.nS_Param{iEq}(i),'%4.4f')]);
            fprintf(1,['\n     Std Error (p. value): ', num2str(Spec_Output.Coeff_SE.nS_Param{iEq}(i),'%4.4f'), ...
                ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.nS_Param{iEq}(i))/Spec_Output.Coeff_SE.nS_Param{iEq}(i),nr-numel(param))),'%4.2f'),')']);
            
        end
    end
    
    if S_Var{iEq}==0
        fprintf(1,'\n\nNon Switching Variance of model ');
        fprintf(1,['\n     Value:                ', num2str(Spec_Output.Coeff.covMat{1}(iEq,iEq),'%4.6f')]);
        fprintf(1,['\n     Std Error (p. value): ', num2str(Spec_Output.Coeff_SE.covMat{1}(iEq,iEq),'%4.4f'), ...
            ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.covMat{1}(iEq,iEq))/Spec_Output.Coeff_SE.covMat{1}(iEq,iEq),nr-numel(param))),'%4.2f'),')']);
    end
    
    switch distrib
        case 't'
            
            if S_df{iEq}==0
                
                fprintf(1,'\nNon Switching Degrees of Freedom (t distribution)');
                fprintf(1,['\n     Value:                ', num2str(Spec_Output.Coeff.df{iEq}(1),'%4.4f')]);
                fprintf(1,['\n     Std Error (p. value): ', num2str(Spec_Output.Coeff_SE.df{iEq}(1),'%4.4f') , ...
                    ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.df{iEq}(1))/Spec_Output.Coeff_SE.df{iEq}(1),nr-numel(param))),'%4.2f'),')']);
            end
            
        case 'GED'
            
            if S_K{iEq}==0
                fprintf(1,'\nNon Switching k parameter (GED distribution)');
                fprintf(1,['\n     Value of k:                 ', num2str(Spec_Output.Coeff.K{iEq}(1),'%4.4f')]);
                fprintf(1,['\n     Std Error (p. value):       ', num2str(Spec_Output.Coeff_SE.K{iEq}(1),'%4.4f') , ...
                    ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.K{iEq}(1))/Spec_Output.Coeff_SE.K{iEq}(1),nr-numel(param))),'%4.2f'),')']);
            end
    end
    
    
    if any([S_Var{iEq},S_df{iEq},S_K{iEq}])
        fprintf(1,'\n\n--->   Switching Parameters (Distribution Parameters)  <---\n');
        
        for j=1:k
            
            fprintf(1,['\nState ', num2str(j)]);
            
            if S_Var{iEq}
                        fprintf(1,['\n    Model''s Variance:      ', num2str(Spec_Output.Coeff.covMat{j}(iEq,iEq),'%4.6f')]);
                        fprintf(1,['\n    Std Error (p. value):  ',num2str(Spec_Output.Coeff_SE.covMat{j}(iEq,iEq),'%4.4f') , ...
                            ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.covMat{j}(iEq,iEq)/Spec_Output.Coeff_SE.covMat{j}(iEq,iEq)),nr-numel(param))),'%4.2f'),')']);
            end
                    
            switch distrib
                    
                case 't'
                    
                    if S_df{iEq}
                        fprintf(1,['\n    Degrees of Freedom (t dist): ', num2str(Spec_Output.Coeff.df{iEq}(1,j),'%4.4f')]);
                        fprintf(1,['\n    Std Error (p. value):        ', num2str(Spec_Output.Coeff_SE.df{iEq}(1,j),'%4.4f') , ...
                            ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.df{iEq}(1,j))/Spec_Output.Coeff_SE.df{iEq}(1,j),nr-numel(param))),'%4.2f'),')']);
                    end
                    
                    
                case 'GED'
                    
                    if S_K{iEq}
                        fprintf(1,['\n    Value of k (GED dist)       ', num2str(Spec_Output.Coeff.K{iEq}(j),'%4.4f')]);
                        fprintf(1,['\n    Std Error (p. value):       ', num2str(Spec_Output.Coeff_SE.K{iEq}(j),'%4.4f') , ...
                            ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.K{iEq}(j))/Spec_Output.Coeff_SE.K{iEq}(j),nr-numel(param))),'%4.2f'),')']);
                    end
                    
            end
        end
    end
    
    fprintf(1,'\n\n--->   Switching Parameters (Regressors)  <---');
    
    if n_S{iEq}==0
        fprintf(1,'\n\nThere was no switching parameters for the regressors in Equation #%i. Skipping this result',iEq);
    else
        
        for i=1:n_S{iEq}
            fprintf(1,['\n\nSwitching Parameters for Equation #' num2str(iEq) ' - Indep column ', num2str(S_S{iEq}(i)),'\n']);
            
            for j=1:k
                fprintf(1,['\nState ', num2str(j)]);
                fprintf(1,['\n   Value:                ', num2str(Spec_Output.Coeff.S_Param{iEq}(i,j),'%4.4f')]);
                fprintf(1,['\n   Std Error (p. value): ', num2str(Spec_Output.Coeff_SE.S_Param{iEq}(i,j),'%4.4f') , ...
                    ' (',num2str(2*(1-tcdf(abs(Spec_Output.Coeff.S_Param{iEq}(i,j))/Spec_Output.Coeff_SE.S_Param{iEq}(i,j),nr-numel(param))),'%4.2f'),')']);
                
            end
        end
    end
end
fprintf(1,'\n\n---> Transition Probabilities Matrix (p-value) <---\n');

% fix for transtion matrix with new algorithm (using fminsearch)
if (strcmp(advOpt.optimizer,'fminsearch'))||(strcmp(advOpt.optimizer,'fminunc'))
    Spec_Output.Coeff_pValues.p=[Spec_Output.Coeff_pValues.p ; repmat(nan,1,k)]; 
end
   
pValue_P=Spec_Output.Coeff_pValues.p;

for i1=1:k
    fprintf(1,'\n      ');
    for i2=1:k
        fprintf(1,'%4.2f (%4.2f)   ',P(i1,i2),pValue_P(i1,i2));
    end
end

% Message warnging about std of p matrix (COMMENTED)
% if strcmp(advOpt.optimizer,'fminsearch')||strcmp(advOpt.optimizer,'fminunc')
%     fprintf(1,['\n* By using fminsearch or fminunc as the estimation function, \n' ...
%         'the standard errors of the last row of the transition matrix  are not available.']);
% end

fprintf(1,'\n\n---> Expected Duration of Regimes <---\n\n');

for i=1:k
    fprintf(1,['     ' 'Expected duration of Regime #%i: %4.2f time periods\n'],i,stateDur(i));
end

if ~advOpt.diagCovMat

    fprintf(1,'\n\n---> Covariance Matrix <---\n');
    for ik=1:k
        fprintf(1,['\nState ', num2str(ik)]);
        pValue_covMat=2*(1-tcdf(abs(Spec_Output.Coeff.covMat{ik})./Spec_Output.Coeff_SE.covMat{ik},nr-numel(param)));

        for iEq=1:nEq
            fprintf(1,'\n      ');
            for jEq=1:nEq
                fprintf(1,'%4.5f (%4.5f,%4.2f)   ',Spec_Output.Coeff.covMat{ik}(iEq,jEq),Spec_Output.Coeff_SE.covMat{ik}(iEq,jEq),pValue_covMat(iEq,jEq));
            end
        end
    end
end

fprintf(1,'\n');