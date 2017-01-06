% Equality and inequality constraints for ms(k)-ar(p) estimation

function [c,ceq]=confuneq_MS_Regress(param)

global k_global; % Getting k
global Coeff_Tag_Global;
global constCoeff_Global;



Coeff_Tag=Coeff_Tag_Global;
k=k_global;
constCoeff=constCoeff_Global;

typeCall='estimation';
[Coeff]=param2spec(param,Coeff_Tag,constCoeff,typeCall);

% Organizing transition probabilities
% 
% for i=0:k-1
%     Coeff.p(i+1,1:k)=param(end-k^2+1+i*k:end-k^2+k+i*k);
% end

% Inequality and equality constrains for the transition probabilities

c=0;
ceq=sum(Coeff.p)-1;