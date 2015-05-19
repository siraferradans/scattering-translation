function [v,meta]=scat2vector(S)

v = [];
%meta = cell(1,size(S,2));
for m=1:size(S,2) 
    %meta{m} = zeros(length(S{m}));
    for r = 1:length(S{m})
        v=[v(:); S{m}{r}.l1(:)];%not taking into account the l2-component!
        
        meta{m}{r}.scale = S{m}{r}.scale;
        meta{m}{r}.orientation = S{m}{r}.orientation;    
    end 
end 

