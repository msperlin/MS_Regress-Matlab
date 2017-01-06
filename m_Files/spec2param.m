% function for converting form

function [Spec_Tag,out_param]=spec2param(Spec)

Fields=fieldnames(Spec);
nField=length(Fields);

count=0;
for n=1:nField
    
    if ~iscell(eval(['Spec.',Fields{n}]))
        str=['Spec.',Fields{n}];
        [nr,nc]=size(eval(str));
        
        for i=1:nr
            for j=1:nc
                count=count+1;
                out_param(count)=eval([str,'(',num2str(i),',',num2str(j),')']);
                
                eval(['Spec_Tag.',Fields{n},'(',num2str(i),',',num2str(j),')','=',num2str(count) ';']);
           
            end
        end
        
    else
        
        str=['Spec.',Fields{n}, '{iCell}' ];
        nCell=numel(eval(['Spec.',Fields{n}]));
        
        for iCell=1:nCell
            
            [nr,nc]=size(eval(str));
            
            if (nr==0)
                eval(['Spec_Tag.',Fields{n}, '{' num2str(iCell) '}','={};']);
            end
            
            for i=1:nr
                for j=1:nc
                    count=count+1;
                    out_param(count)=eval([str,'(',num2str(i),',',num2str(j),')']);
                    eval(['Spec_Tag.',Fields{n},'{' num2str(iCell) '}','(',num2str(i),',',num2str(j),')','=',num2str(count) ';']);
                    
                end
            end
        end
    end
end

