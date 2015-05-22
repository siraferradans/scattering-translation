%% Testing the learning of a translation matrix

%First: generate a DB 
%clear all
J1=8;
Q=8;
J2=8;
L=8;
num_signals=(J1*Q+1);
% m=2*(J1*Q+1);
% n=4*(J2*L+1);
 
[x,y,meta,filt2d,opt2d,filt1d,opt1d]=generate_Sound2Image_DB(J1,Q,J2,L,num_signals);

T = linsolve(x,y);

save('./Synth_DB_small.mat','T','x','y','meta','filt2d','opt2d','filt1d','opt1d');

% 
%% Testing
opt1d.l2scatt=0;
playsound=@(x)soundsc(fftshift(x-mean(x(:))),1024);
audio=[];
for indx = 1:Q:size(x,1)
    
    Sa=vector2scat(x(indx,:),meta,[2 1]);
    [audio(end+1,:),~]= newscatt_synthesis_mgrid(Sa, filt1d, opt1d, x(indx,:), max(x(indx,:)));
   
    yy = T'*x(indx,:)';

    %Sa=vector2scat(x(1,:),meta,[2 1]);

%    %    to compare with the original
%     Si=vector2scat(y(indx,:),meta,[2 2]);
%     opt2d.l2scatt=false;
%     [reco1,energy1]= newscatt_synthesis_mgrid(Si, filt2d, opt2d, x(indx,:), max(x(indx,:)));

    Sii=vector2scat(yy',meta,[2 2]);
    opt2d.l2scatt=false;
    [reco2,energy2]= newscatt_synthesis_mgrid(Sii, filt2d, opt2d, x(indx,:), max(x(indx,:)));

     figure;imshow(reco2,[])
     playsound(audio(end,:))
     hold on;
end



%%%%%%%%%%%%%%%%%%%%
%% testing with given parameters
f = ones(2^J1,1);
j=2;q=8;
signal = real(ifft(f.*filt1d{1}.psi{1}{(j-1)*Q+q}{1}));



plot(fftshift(signal));
fs = 22050;
%soundsc(fftshift(signal),fs)
        %%
bb = zeros(22050,1);
bb=signal(1:22050,1);
audiowrite('./target.wav',target,22050);

% target = signal-min(signal(:));
% target = target/norm(target(:));

Sa= fwdscatt(signal,filt1d, opt1d);
a = scat2vector(Sa);
i = T'*a;
Si=vector2scat(i,meta,[2 2]);

opt2d.l2scatt = 0;
[reco1,energy1]= newscatt_synthesis_mgrid(Si, filt2d, opt2d, signal, max(signal(:)));
figure;imshow(reco1,[])
