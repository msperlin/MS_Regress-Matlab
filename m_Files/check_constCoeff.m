if any([~isfield(advOpt.constCoeff,'nS_Param'), ...
        ~isfield(advOpt.constCoeff,'S_Param') , ...
        ~isfield(advOpt.constCoeff,'covMat')  , ...
        ~isfield(advOpt.constCoeff,'p') ])
    
    str=sprintf(['In the construction of constCoeff, its missing one (or more) of the fields:\n' ...
        'nS_param\n','S_Param\n','Std\n','p\n','See Example files (and pdf) for details of how to use advOpt.constCoeff']);
    error(str);
end

if any([~iscell(advOpt.constCoeff.nS_Param) , ...
        ~iscell(advOpt.constCoeff.S_Param)  , ...
        ~iscell(advOpt.constCoeff.covMat)   , ...
        ~iscell(advOpt.constCoeff.p)])
    
    str=sprintf(['In the construction of constCoeff, all members should be cell arrays ' ...
        'with number of elements equal to the number of equations in system (see example files']);
    error(str);
    
end

switch distrib
    case 't'
        if ~isfield(advOpt.constCoeff,'df')
            error('In argument advOpt.constCoeff, its missing the parameter df for the t distribution')
        end
        
    case 'GED'
        if ~isfield(advOpt.constCoeff,'K')
            error('In argument advOpt.constCoeff, its missing the parameter K for the GED distribution')
        end
end

%  testing for probabilities in transition matrix

for jk=1:k
    for ik=1:k
        myTest(ik,jk)=strcmp(advOpt.constCoeff.p{ik,jk},'e');
    end
    
    if numel(unique(myTest(:,jk)))>1
        error(['In argument advOpt.constCoeff.p, each column must have either probabitities ' ...
            'summing to 1 or ''e'' strings. This means that for transition matrix you cannot ' ...
            'have mixed "estimate this" flag (''e'') and probabilities in same state']) ;
    end
end

% Making sure that probabilities in transition matrix sum to 1

if any(any(myTest==0))
    myIdx=find(myTest(1,:)==0);
    
    for i=1:numel(myIdx)
        for ik=1:k
            tempPmatrix(ik,myIdx)=advOpt.constCoeff.p{ik,i};
        end
    end
    
    if any( (sum(tempPmatrix)>1.0001)|(sum(tempPmatrix)<.9999) )
        error(['When setting arbitrary transition probabilities the sum of each collum ' ...
            'in advOpt.constCoeff.p should be equal to 1 (they represent probabilities of a full process.)']);
    end
    
    if any([tempPmatrix<0 tempPmatrix>1])
        error('Any value set in advOpt.constCoeff.p should be lower than 1 and higher than 0.');
    end
end

if sum(sum(strcmp('e',advOpt.constCoeff.p)))~=k*k
   
    if (~strcmp(advOpt.optimizer,'fmincon'))
        str='The use of numerical restrictions in the transition matrix is only possible when using advOpt.optimizer=''fmincon''.';
        error(str);    
    end
end
