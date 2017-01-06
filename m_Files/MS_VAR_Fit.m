% Function for estimation of a Markov Switching Vector Autoregressive Model
%
%   Input:  dep     - Dependent Variable (matrix)
%           
%           nLag    - number of lags in the system
%           k       - Number of States (integer higher than 2)
%           intercept - add intercept ? (1 or 0)
%           advOpt  - A structure with advanced options for algorithm.
%                     See pdf file for details.
%
%   Output: Spec_Output - A structure with all information regarding the
%                         model estimated from the data (see pdf for details).
%
%   Author: Marcelo Perlin (Reading University-ICMA/UK)
%   Contact:  m.perlin@icmacentre.ac.uk

function [Spec_Output]=MS_VAR_Fit(dep,nLag,k,intercept,advOpt)

if nargin()==3
    intercept=1;
    advOpt.distrib='Normal';
    advOpt.std_method=1;
    advOpt.optimizer='fminsearch';
end

if nargin()==4
    advOpt.distrib='Normal';
    advOpt.std_method=1;
    advOpt.optimizer='fminsearch';
end

if ~isfield(advOpt,'optimizer')
    advOpt.optimizer='fminsearch';
end

if nLag<0
    error('Input nLag should be a positive integer')
end

nEq=size(dep,2);
nRow=size(dep,1);

% building lagged variables 

if intercept
    indep=ones(nRow,1);
else 
    indep=[];
end
    
for iEq=1:nEq
    for j=1:nLag
        temp=[zeros(j,1) ; dep(1:end-j,iEq)];
        indep=[indep,temp];
    end
end

switch advOpt.distrib
    case 'Normal'
        n_dist_param=1;
    case 't'
        n_dist_param=2; % never used (MV only works for NORMAL)
    case 'GED'
        n_dist_param=2;  % never used (MV only works for NORMAL)
end

for iEq=1:nEq
    S{iEq}=repmat(1,1,size(indep,2)+n_dist_param);
end

advOpt.printOut=0;
[Spec_Output]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

doOutScreen_MSVAR();

