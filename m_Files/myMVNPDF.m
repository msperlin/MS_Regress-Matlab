function pdf=myMVNPDF(dep,Cond_mean,covMat)

% [n,d]=size(dep); 
[nr,nc]=size(dep); 

% S=covMat;
% X=dep;
% m=Cond_mean;
% 
% den = (2*pi*det(S))^(d/2); 
% Xn = X-m(ones(n,1),:); 
% pdf = zeros(n,1); 
%   
%  % new and fast call 
% pdf = exp(-0.5*sum((Xn(:,:)/S).*Xn(:,:) ,2))/den; 
% n=pdf;

pdf=1/(((2*pi())^(nc/2))*sqrt(det(covMat))).* ...
        exp(-0.5.*sum((dep(:,:)-Cond_mean(:,:))*inv(covMat).*(dep(:,:)-Cond_mean(:,:)),2));

% for i=1:nr
%     pdf(i,1)=1/(((2*pi())^(nc/2))*sqrt(det(covMat))).* ...
%         exp(-0.5.*(dep(i,:)-Cond_mean(i,:))*(inv(covMat)*(dep(i,:)-Cond_mean(i,:))'));
% end

% S1 = inv(S); 
% for ix=1:n 
%    pdf2(ix) = exp(-0.5*Xn(ix,:)*S1*(Xn(ix,:).'))/den; 
% end 
% 
% for i=1:n 
%    pdf3(i) = exp(-0.5*Xn(i,:)*S1*(Xn(i,:).'))/den; 
% end 
% [nr nc]=size(dep);
% 
% den=(2*pi*det(covMat))^(nc/2);
% 
% n = exp(-0.5*sum((dep(:,:)/covMat).*dep(:,:) ,2))/den; 

% for i=1:nr
%     n(i,1)=(exp(-1/2.*[dep(i,:)-Cond_mean(i,:)]'*inv(covMat)*[dep(i,:)-Cond_mean(i,:)]))/den;
% end
