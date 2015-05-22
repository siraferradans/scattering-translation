function [x,y,meta,filt2d,opt2d,filt1d,opt1d]=generate_Sound2Image_DB(J1,Q,J2,L,num_signals)
%num_signals: number of instances we want in the data_set
%J1: 2^J1 size of the audio signal
%Q: number of wavelets per octave
%J2: (2^J2,2^J2) size of the image signal
%L: number of angles to analyze 
%
% Output: 
% x: size (num_signals,J1 Q)
% y: size (num_signals,J2 L)

%% start by the image part
% Compute the Scattering of every IMAGE and save it
Ni = 2^J2;
Ji = {J2,1,1};
Li = {L,1,1};
Qi = {1,1,1};
onedim=false;
[filt2d,opt2d]=getTransfilters_for_reconstruction(Ni,Ji,Qi,Li,onedim);

%% Audio part
Na = 2^J1;
Ja = {J1,1,1};
La = {1,1,1};
Qa = {Q,1,1};
onedim=true;
[filt1d,opt1d]=getTransfilters_for_reconstruction(Na,Ja,Qa,La,onedim);


%% generate the scatt and save the pairing
signal_base = rand(Na,1);
Q = opt1d.Q1;

%% pairs of j,q
%get all possible subsets of Q and J
%qsubset = logical(dec2bin(0:2^Q-1)-'0');
%qsubset = qsubset(2:end,:);%take out the 'all 0s'
%qsubset = qsubset(randperm((2^Q-1)),:); %in case we stop before all the possible perm. are done.

qsubset = eye(Q,Q);
%qsubset = cat(1,eye(Q,Q),qsubset);

%jsubset = logical(dec2bin(0:2^(J1)-1)-'0');
%jsubset = jsubset(2:end,:);%take out the 'all 0s'
%jsubset = jsubset(randperm((2^J1-1)),:);
%jsubset = cat(1,eye(J1,J1),jsubset);

jsubset = eye(J1,J1);

x = [];%zeros(0,(1+J1*Q)*2);  
y = [];%zeros(0,(1+J2*L)*4);

k=1;
for codej = 1:size(jsubset,1) %for all subsets of J
    js= find(jsubset(codej,:));
    for codeq=1:size(qsubset,1) %for all subsets of Q
        qs=find(qsubset(codeq,:));
        
        xx = zeros(size(signal_base));
        for j=js
            for q = qs
                f = filt1d{1}.psi{1}{(j-1)*Q+q}{1};
                xx = xx+(signal_base).*f;
            end
        end
        
        target = real(ifft(xx));
        
        target = target-min(target(:));
        target = target/norm(target);
        
        Sa= fwdscatt(target,filt1d, opt1d);
        Si= Scataudio_to_Scatimage(Sa,filt1d,opt2d.J1);
        
        %     [reco,energy]= newscatt_synthesis_mgrid(Si, filt_2d, opt_2d, target, max(target(:)));
        %     %see difference
        %     [Sr,~]= fwdscatt(reco,filt_2d, opt_2d);%for debug
        %
        y(end+1,:) = scat2vector(Si);
        x(end+1,:) = scat2vector(Sa);
        
        
        c = hsv2rgb(rand(1,3));
        subplot(1,2,1);plot(y(k,:),'Color',c);hold on;
        subplot(1,2,2);plot(x(k,:),'Color',c);hold on;
        k = k+1;
        if k>num_signals
            [~,meta]= scat2vector(Si);
            return;
        end
    end
end

[~,meta]= scat2vector(Si);