
%% Test Gaussian
%clear all; 
close all

% select=@(target,r)target(1:end/r);
% path='../../data/mcdermott/';
% loadresample=@(file,r)select(resample(audioread([path file]),1,r),r);

% fl=@(x)x(:);
% remap=@(x,m,M)(M-m)*(x-min(fl(x)))/(max(fl(x))-min(fl(x)))+m;
% remap01=@(x)remap(x,0,1);

N = 256;
% target = loadresample('orig_ex_1.wav',2);
% target = target(1:N);
%for audio!!

options.N = 256;
options.localized=1;
options.J1=min(log2(options.N),12);
options.J2=min(log2(options.N),12);
options.J3=min(log2(options.N),12);
options.L1=1;
options.L2=1;
options.L3=1;
options.Q1=8;
options.niters=1000;
options.dataset='unidim';
options.lambda=+1e-6;
options.multigrid=0;
options.recenter = 0;
options.lr0=0.1;%8e-1;
options.init_with_first=0;
options.init_with_sec=0;
options.maxorder=1; % !!!!
options.onedim=1;
options.os=2;

mmax = getoptions(options,'maxorder',2);
init_with_first = getoptions(options,'init_with_first',0);%set to 1 to initialize the full reconstruction with first order reconstruction.
init_with_sec = getoptions(options,'init_with_sec',0);%set to 1 to initialize the full reconstruction with first order reconstruction.
fouriereq=getoptions(options,'fouriereq',0);%replace scattering with spectral equalization
N  = getoptions(options,'N',1024);
J1 = getoptions(options,'J1',8);
dataset=getoptions(options,'dataset','unidim');
Q1 = getoptions(options,'Q1',1);
J2 = getoptions(options,'J2',8);
Q2 = getoptions(options,'Q2',1);
J3 = getoptions(options,'J3',8);
onlyfirst = getoptions(options,'onlyfirst',0);%set to 1 for only first order scattering
localized = getoptions(options,'localized',1);%set to 1 to output localized scattering coefficients (at scale 2^J1).
multigrid = getoptions(options,'multigrid',1);%set to 1 for multigrid reconstruction
lr0 = getoptions(options,'lr0',1e-1);%learning rate
recenter = getoptions(options,'recenter',0);%set to 1 to recenter reconstructions in the case of delocalized scattering
l2scatt = getoptions(options,'l2scatt',0);%set to 1 to include also l2 norms along with l1 norms
niters = getoptions(options,'niters',200);
fouriereq=getoptions(options,'fouriereq',0);%replace scattering with spectral equalization

onedim = getoptions(options,'onedim',1);

options.lr=lr0;%2e-2;
options.momentum=0.9;
options.niters=niters;
options.rhotol=0.01;
options.mirror=0;%do not use mirror symmetry
options.periodinput=0;
options.border=0;%do not remove borders
options.multigrid=multigrid;
options.usepinv=1;% use dual filters in the reconstruction step.
options.N=N;
options.J1=J1;
options.Q1=Q1;
options.J2=J2;
options.Q2=Q2;
options.J3=J3;
options.splines=0;


options.l2scatt=1;
options.positive = 0;
options

%% Filter generation
% audio part
[filters,lpal] = generate_scatt_filters(options);
J = log2(N);
%number of coefficients:
1+options.Q1*options.J1
%1+options.Q1*options.J1+(options.Q1*options.Q2*options.J1*(options.J2-1))/2

% image part
opt_2d = options;
opt_2d.L1=8;
opt_2d.L2=1;
opt_2d.Q1=1;
opt_2d.Q2=1;

opt_2d.onedim=false;
opt_2d.softthreshold=0;
opt_2d
[filt_2d,lpal_2d] = generate_scatt_filters(opt_2d);

1+opt_2d.L1*opt_2d.J1%+(opt_2d.L1*opt_2d.J1*opt_2d.L2*(opt_2d.J2-1))/2

%%%%%%%%%%%%%%%% Application to signals %%%%%%%%%%%%%%%%%%%

%% Gaussian
target = randn(N,1);
target = target-min(target(:));
target = target/norm(target);

Sa= fwdscatt(target,filters, options);

%%%%%% testing...
%Ojo Q=1!! sino no funciona!  y los valores no pueden ser negativos 
% reco= newscatt_synthesis_mgrid(Sa, filters, options, target, max(target(:)));
% a = scat2vector(Sa);
% [Saux,~]= fwdscatt(reco,filters, options);
% r = scat2vector(Saux);
% plot(a);hold on; plot(r,'r')
%%%%%

% translate Sa into a Simages
% just taking into account the first and second layer, but we need to
% change the meta data

Si=Scataudio_to_Scatimage(Sa,filters,opt_2d.J1);
[reco,energy]= newscatt_synthesis_mgrid(Si, filt_2d, opt_2d, target, max(target(:)));

%see difference
[Sr,~]= fwdscatt(reco,filt_2d, opt_2d);%for debug
figure;
r = scat2vector(Sr);
s = scat2vector(Si);
plot(r);hold on; plot(s,'r')
figure;subplot(1,2,1);hist(reco(:),512);title(['M=' num2str(mean(reco(:))) ' std=' num2str(std(reco(:)))])
subplot(1,2,2);hist(target(:),512);title(['M=' num2str(mean(target(:))) ' std=' num2str(std(target(:)))])

%% Poisson
target = zeros(N,1);
target([20 56 100 112 200]) = 1;%poissrnd(0.5mmm,N,1);
target = target-min(target(:));
target = target/norm(target);

Sa= fwdscatt(target,filters, options);
Si=Scataudio_to_Scatimage(Sa,filters,opt_2d.J1);
[reco,energy]= newscatt_synthesis_mgrid(Si, filt_2d, opt_2d, target, max(target(:)));

%see difference
[Sr,~]= fwdscatt(reco,filt_2d, opt_2d);%for debug
figure;
r = scat2vector(Sr);
s = scat2vector(Si);
plot(r);hold on; plot(s,'r')
figure;subplot(1,2,1);hist(reco(:),512);title(['M=' num2str(mean(reco(:))) ' std=' num2str(std(reco(:)))])
subplot(1,2,2);hist(target(:),512);title(['M=' num2str(mean(target(:))) ' std=' num2str(std(target(:)))])


% %% Synthesizing new Gaussian texture 1d
% cov=@(f) abs( fft(f-mean(f(:))).^2 ) ;
% synthesisGaussian=@(N,C,m)1/sqrt(options.N)*real(ifft(sqrt(C).*fft(randn(N,1))))+m;
% reco_Gaussian = synthesisGaussian(options.N,cov(target),mean(target));
% 
% subplot(131);plot(target);title('Original')
% subplot(132);plot(reco);title('Reco. Scattering')
% subplot(133);plot(reco_Gaussian);title('Gaussian')

