nr=size(dep,1);
param=Spec_Output.param;

fprintf(1,'\n\n***** Numerical Optimization Converged *****\n\n');
fprintf(1,['Final log Likelihood: ',num2str(Spec_Output.LL),'\n']);
fprintf(1,['Number of estimated parameters: ',num2str(Spec_Output.Number_Parameters),'\n']);
fprintf(1,['Number of Observations: ',num2str(Spec_Output.nObs),'\n']);
fprintf(1,['Number of Equations: ',num2str(Spec_Output.nEq),'\n']);
fprintf(1,['Optimizer: ',Spec_Output.advOpt.optimizer,'\n']);
fprintf(1,['Number of Equations in System: ',num2str(nEq),'\n']);
fprintf(1,['Distribution Assumption -> ',Spec_Output.advOpt.distrib,'\n']);
fprintf(1,['Standard error calculation -> ',num2str(Spec_Output.advOpt.std_method) '\n']);

for iEq=1:nEq
    fprintf(1,['\n***** Final Parameters for Equation #' num2str(iEq) ' ***** \n\n']);
    
    if intercept    % display intercept
        fprintf(1,[blanks(5) 'Intercept - Parameter Value (Standard Error, p value)\n']);
        
        for ik=1:k
            seValue=Spec_Output.Coeff_SE.S_Param{iEq}(1,ik);
            pValue=2*(1-tcdf(abs(Spec_Output.Coeff.S_Param{iEq}(1,ik))/seValue,size(dep,1)-numel(Spec_Output.param)));
            fprintf(1,[blanks(10) 'State %i, Intercept = %4.2f (%4.2f,%4.2f) \n'],ik,Spec_Output.Coeff.S_Param{iEq}(1,ik),seValue,pValue);
        end
    end
    
    myCounter=intercept+1;
    for jEq=1:nEq
        fprintf(1,[blanks(5) 'Dependent Variable #%i - Parameter Value (Standard Error, p value)\n'],jEq);
        
        for iLag=1:nLag
            for ik=1:k
                seValue=Spec_Output.Coeff_SE.S_Param{iEq}(myCounter,ik);
                pValue=2*(1-tcdf(abs(Spec_Output.Coeff.S_Param{iEq}(myCounter,ik))/seValue,size(dep,1)-numel(Spec_Output.param)));
                fprintf(1,[blanks(10) 'State %i, Lag %i = %4.2f (%4.2f,%4.2f) \n'],ik,iLag,Spec_Output.Coeff.S_Param{iEq}(myCounter,ik),seValue,pValue);
                
            end
            myCounter=myCounter+1;
        end
    end
    
end
fprintf(1,'\n---> Transition Probabilities Matrix (p-value) <---\n');

% fix for transtion matrix with new algorithm (using fminsearch)
if (strcmp(advOpt.optimizer,'fminsearch'))||(strcmp(advOpt.optimizer,'fminunc'))
    Spec_Output.Coeff_pValues.p=[Spec_Output.Coeff_pValues.p ; repmat(nan,1,k)]; 
end
   
pValue_P=Spec_Output.Coeff_pValues.p;

for i1=1:k
    fprintf(1,'\n      ');
    for i2=1:k
        fprintf(1,'%4.2f (%4.2f)   ',Spec_Output.Coeff.p(i1,i2),pValue_P(i1,i2));
    end
end

fprintf(1,'\n\n---> Expected Duration of Regimes <---\n\n');

for i=1:k
    fprintf(1,['     ' 'Expected duration of Regime #%i: %4.2f time periods\n'],i,Spec_Output.stateDur(i));
end

fprintf(1,'\n---> Covariance Matrix <---\n');
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

disp(' ');