function [v,meta]=scat2vector(S)

v = [];
for m=1:size(S,2) 
    for r = 1:length(S{m})
        v=[v(:); S{m}{r}.l1(:)];%; S{m}{r}.l2(:)];
        
        meta{m}{r}.scale = S{m}{r}.scale;
        meta{m}{r}.orientation = S{m}{r}.orientation;    
    end 
end 

