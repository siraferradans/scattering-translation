function S=vector2scat(v,meta,dim)

len = dim(1)*dim(2);

k = 0;
for m=1:size(meta,2) 
    
    for r = 1:length(meta{m})
        
        S{m}{r}.l1 = reshape(v( len*k+1:len*(k+1) ),dim(1),dim(2));
        k=k+1;
%         S{m}{r}.l2 = reshape(v( len*k+1:len*(k+1) ),dim(1),dim(2));
%         k=k+1;
%         
        S{m}{r}.scale = meta{m}{r}.scale;
        S{m}{r}.orientation = meta{m}{r}.orientation;    
    end 
end 
