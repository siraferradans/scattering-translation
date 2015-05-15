close all;
clear all;

options.N = 256;
options.localized=1;
options.J1=min(log2(options.N),12);
options.J2=min(log2(options.N),12);
options.J3=min(log2(options.N),12);
options.L1=1;
options.L2=1;
options.L3=1;
options.Q1=1;
options.niters=1000;
options.dataset='unidim';
options.lambda=+1e-6;
options.multigrid=1;
options.recenter = 0;
options.lr0=0.1;%8e-1;
options.init_with_first=0;
options.init_with_sec=0;
options.maxorder=1;
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
options.splines=1;


options.l2scatt=1;
options.multigrid=1;
options.positive = 0;
options

if options.positive
options.softthreshold=1;
removesign=@(x)(abs(x));
else
options.softthreshold=0;
removesign=@(x)(x);
end

[filters,lpal] = generate_scatt_filters(options);


target=rand(N,1);%get signal as instance
target = removesign(target)/norm(target);

Sa= fwdscatt(target,filters, options);
% 
% init = removesign(randn(N,1));
% init = init/norm(init);
% options.input = init;

[reco, E]= newscatt_synthesis_mgrid(Sa, filters, options, target, max(abs(target(:))));

e = scat2vector(Sa);
Sr= fwdscatt(reco,filters, options);
r = scat2vector(Sr);

plot(e);hold on;plot(r,'r')



%%
opt_2d = options;
opt_2d.L1=1;
opt_2d.L2=1;
opt_2d.Q1=1;
opt_2d.Q2=1;

opt_2d.onedim=false;

[filt_2d,lpal_2d] = generate_scatt_filters(opt_2d);

1+opt_2d.L1*opt_2d.J1

opt_2d
% translate Sa into a Simages
% just taking into account the first and second layer, but we need to
% change the meta data

Si=Scataudio_to_Scatimage(Sa,filters,opt_2d.J1);

reco= newscatt_synthesis_mgrid(Si, filt_2d, opt_2d, target, max(target(:)));
[Sr,~]= fwdscatt(reco,filt_2d, opt_2d);%for debug
