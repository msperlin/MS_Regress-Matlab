function [Spec_Out]=param2spec(param,Spec_Tag,constCoeff,typeCall)

Fields=fieldnames(Spec_Tag);
nField=length(Fields);
nEq=size(Spec_Tag.covMat{1},1);

count=0;
for n=1:nField
    
    if ~iscell(eval(['Spec_Tag.',Fields{n}]))
        
        str=['Spec_Tag.',Fields{n}];
        
        [nr,nc]=size(eval(str));
        
        for i=1:nr
            for j=1:nc
                count=count+1;
                
                spec_chunk=['Spec_Tag.',Fields{n},'(',num2str(i),',',num2str(j),')'];
                str1=['Spec_Out.',Fields{n},'(',num2str(i),',',num2str(j),')','='];
                
                
                str2=['param(',spec_chunk,');'];
                
                if eval([spec_chunk '==0;'])
                    switch typeCall
                        case 'estimation'
                            eval([str1 'constCoeff.' Fields{n} '{',num2str(i),',',num2str(j),'};']);
                        case 'se_calculation'
                            eval([str1 'NaN;']);
                    end
                    
                    continue;
                end
                
                eval([str1 str2]);
                
                %                 out_param(count)=eval([str,'(',num2str(i),',',num2str(j),')']);
                %                 eval(['Spec_Tag.',Fields{n},'(',num2str(i),',',num2str(j),')','=',num2str(count)]);
            end
        end
        
        
    else
        
        nCell=numel(eval(['Spec_Tag.',Fields{n}]));
        
        str=['Spec_Tag.',Fields{n},'{iCell}'];
        for iCell=1:nCell
            
            [nr,nc]=size(eval(str));
            
            for i=1:nr
                for j=1:nc
                    count=count+1;
                    
                    
                    spec_chunk=['Spec_Tag.',Fields{n},'{iCell}','(',num2str(i),',',num2str(j),')'];
                    str1=['Spec_Out.',Fields{n},'{iCell}','(',num2str(i),',',num2str(j),')','='];
                    
                    str2=['param(',spec_chunk,');'];
                    
                    if eval([spec_chunk '==0;'])
                        switch typeCall
                            case 'estimation'
                                eval([str1 'constCoeff.' Fields{n} '{iCell}' '{',num2str(i),',',num2str(j),'};']);
                            case 'se_calculation'
                                eval([str1 'NaN;']);
                        end
                        
                        continue;
                    end
                    
                    eval([str1 str2]);
                    
                    %                 out_param(count)=eval([str,'(',num2str(i),',',num2str(j),')']);
                    %                 eval(['Spec_Tag.',Fields{n},'(',num2str(i),',',num2str(j),')','=',num2str(count)]);
                end
            end
            
        end
    end
end
