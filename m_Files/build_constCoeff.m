% building constCoeff

if nargin()==4||(isfield(advOpt,'constCoeff')==0) % if number of arg < 4, build constraint structure
                                                  % and fill it with 'e'
    
    for ik=1:k
        advOpt.constCoeff.covMat{1,ik}=cell(nEq);
    end
    
    advOpt.constCoeff.p=cell(k,k);
    
    for ik=1:k
        for i=1:nEq
            for j=1:nEq
                advOpt.constCoeff.covMat{ik}{i,j}='e';
            end
        end
    end
    
    for iEq=1:nEq
        
        switch distrib
            case 't'
                advOpt.constCoeff.df{iEq}=cell(1,S_df{iEq}*(k)+(1-S_df{iEq}));
                
                for i=1:size(advOpt.constCoeff.df{iEq},2)
                    advOpt.constCoeff.df{iEq}{1,i}='e';
                end
                
            case 'GED'
                advOpt.constCoeff.K{iEq}=cell(1,S_K{iEq}*(k)+(1-S_K{iEq}));
                
                for i=1:size(advOpt.constCoeff.K{iEq},2)
                    advOpt.constCoeff.K{iEq}{1,i}='e';
                end
        end
        
        advOpt.constCoeff.nS_Param{iEq}=cell(sum(S{iEq}(1:end-n_dist_param)==0),1);
        advOpt.constCoeff.S_Param{iEq}=cell(sum(S{iEq}(1:end-n_dist_param)==1),k*(any(sum(S{iEq}(1:end-n_dist_param)==1))));
        
        for i=1:size(advOpt.constCoeff.nS_Param{iEq},1)
            advOpt.constCoeff.nS_Param{iEq}{i}='e';
        end
        
        for i=1:size(advOpt.constCoeff.S_Param{iEq},1)
            for j=1:size(advOpt.constCoeff.S_Param{iEq},2)
                advOpt.constCoeff.S_Param{iEq}{i,j}='e';
            end
        end
        
        for i=1:k
            for j=1:k
                advOpt.constCoeff.p{i,j}='e';
            end
        end
        
        if all(S{iEq}(1:end-n_dist_param)==1)
            advOpt.constCoeff.nS_Param{iEq}{1}=0;
        end
        
        if all(S{iEq}(1:end-n_dist_param)==0)
            advOpt.constCoeff.S_Param{iEq}{1}=0;
        end
    end
end

if advOpt.diagCovMat==1
    for ik=1:k
        for i=1:nEq
            for j=1:nEq
                if i~=j
                    advOpt.constCoeff.covMat{ik}{i,j}=0;
                end
            end
        end
    end
else
    for ik=1:k
        for i=1:nEq
            for j=1:nEq
                if i<j
                    advOpt.constCoeff.covMat{ik}{i,j}=NaN;
                end
            end
        end
    end
end

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
                advOpt.constCoeff.covMat{ik}{iEq,iEq}=NaN; % for cases where variance doesnt switch states
            end
    end
end

switch distrib % not really applicable for iEq iteration (nEq for t is always =1) but keep for future reference
    case 't'
        for ik=2:k
            for iEq=1:nEq
                if S{iEq}(end)==0
                    advOpt.constCoeff.df{iEq}{1,ik}=NaN; % for cases where variance doesnt switch states
                end
            end
        end

    case 'GED'
        for ik=2:k
            for iEq=1:nEq
                if S{iEq}(end)==0
                    advOpt.constCoeff.K{iEq}{1,ik}=NaN; % for cases where variance doesnt switch states
                end
            end
        end
end