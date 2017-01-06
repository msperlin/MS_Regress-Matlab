% checkSize_constCoeff

constCoeff=advOpt.constCoeff;

for i=1:k
    if any(size(constCoeff.covMat{i})~=[nEq nEq])
        error('The argument constcoeff.covMat{i} should have size [k,k].');
    end
end

if any(size(constCoeff.p)~=[k k])
    error('The argument constCoeff.p should be a 2d matrix with size [k k]')
end

for iEq=1:nEq

    if size(constCoeff.nS_Param{iEq},2)>1
        error('The argument constcoeff.nS_Param should be a row vector, not a column vector (or matrix)');
    end


    if ~all(S{iEq}(1:end-n_dist_param)==1)
        if size(constCoeff.nS_Param{iEq},1)~=sum(S{iEq}(1:end-n_dist_param)==0)
            error('The argument constcoeff.nS_Param{iEq} should be a vector with number of rows equal to the number of non switching parameters (indep matrix)')
        end

        for i=1:size(constCoeff.nS_Param{iEq},1)
            if ~isnumeric(constCoeff.nS_Param{iEq}{i})
                if  ~(strcmp(constCoeff.nS_Param{iEq}{i},'e'))
                    error(['Error at constCoeff.nS_Param{iEq}{' num2str(i) '}. Such input can only take a numeric value or string ''e''. See example files for details.']);
                end
            end
        end
    end

    if ~all(S{iEq}(1:end-n_dist_param)==0)
        if size(constCoeff.S_Param{iEq},1)~=sum(S{iEq}(1:end-n_dist_param)==1)
            error('The argument constcoeff.S_Param should have number of rows equal to the number of switching parameters (indep matrix)');
        end

        if size(constCoeff.S_Param{iEq},2)~=(any(S{iEq}(1:end-n_dist_param)==1))*k
            error('The argument constcoeff.S_Param should have number of columns equal to the number of states (k)');
        end
    end

    % checking if fields contain numeric or 'e', only

    for i=1:size(constCoeff.S_Param{iEq},1)
        for j=1:size(constCoeff.S_Param{iEq},2)
            if ~isnumeric(constCoeff.S_Param{iEq}{i,j})
                if  ~(strcmp(constCoeff.S_Param{iEq}{i,j},'e'))
                    error(['Error at constCoeff.S_Param{iEq}{' num2str(i) ',' num2str(j) '}. Such input can only take a numeric value or string ''e''. See example files for details.']);
                end
            end
        end
    end
end

for i=1:k
    for j=1:k
        if ~isnumeric(constCoeff.p{i,j})
            if  ~(strcmp(constCoeff.p{i,j},'e'))
                error(['Error at constCoeff.p{' num2str(i) ',' num2str(j) '}. Such input can only take a numeric value or string ''e''. See example files for details.']);
            end
        end
    end
end
