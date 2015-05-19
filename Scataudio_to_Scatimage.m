function Si=Scataudio_to_Scatimage(Sa,filt_a,J)

Q = filt_a{1}.Q;
%we want to map q->l and j->j
for m=1:size(Sa,2) %m==2 means we only have 1 level {0},{1}
    for r = 1:size(Sa{m},2)
        Si{m}{r}.l1 = cat(2,Sa{m}{r}.l1,Sa{m}{r}.l1);
        Si{m}{r}.l2 = cat(2,Sa{m}{r}.l2,Sa{m}{r}.l2);

        [j,q] = ind2sub([J Q], Sa{m}{r}.scale(m));
        
        if m==2 
            Si{m}{r}.scale = [0 j];
            Si{m}{r}.orientation = [0 q];	
        else
            Si{m}{r}.scale = j;
            Si{m}{r}.orientation = 0;	
        end    
    end 
end 

