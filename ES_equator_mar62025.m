
%% Dan Whitt (daniel.b.whitt@nasa.gov) Jul 11 2024
%% Copyright Dan Whitt
%% Written With Matlab v2023b Fig 6 and 15
%% depednencies:
% ~/OneDrive - NASA/documents/mitgcm_grid.mat
% ~/OneDrive - NASA/documents/mitgcm_RF.mat
% ~/OneDrive - NASA/documents/mitgcm_grid_3.mat
% ~/OneDrive - NASA/documents/mitgcm_grid_2.mat
% ~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc
% ~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc
% ~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc
% ~/OneDrive - NASA/documents/mitgcm_offline_salt_flux_avg_anoms.mat
% ~/OneDrive - NASA/documents/mitgcm_suf_avg.nc
% ~/OneDrive - NASA/documents/mitgcm_hb_avg.nc
clear all
%close all
restoredefaultpath;
addpath ~/'OneDrive - NASA'/software/rainevent-master/seawater/
addpath ~/'OneDrive - NASA'/software/rainevent-master/cmocean/

% set up mit grid
load('~/OneDrive - NASA/documents/mitgcm_grid.mat','DRF','RAC','hFacC');
load('~/OneDrive - NASA/documents/mitgcm_RF.mat','RF');
load('~/OneDrive - NASA/documents/mitgcm_grid_3.mat','DXC','DXG','DYC','DYG');
load('~/OneDrive - NASA/documents/mitgcm_grid_2.mat','RAS','RAW','hFacS','hFacW');

RAC=RAC(41:1460,41:440);
hFacC=hFacC(41:1460,41:440,1:136);
RF=repmat(RF,[1420 400 1]);
RF=RF(:,:,1:136);
DRF=repmat(DRF,[1420 400 1]);
DRF=DRF(:,:,1:136);
RAS=RAS(41:1460,41:440);
RAW=RAW(41:1460,41:440);
hFacW=hFacW(41:1460,41:440,1:136);
hFacS=hFacS(41:1460,41:440,1:136);
DXC=double(DXC(41:1460,41:440));
DXG=DXG(41:1460,41:440);
DYC=DYC(41:1460,41:440);
DYG=DYG(41:1460,41:440);


% constants
g = 9.81; %gravity
rhoref = 1035; % rho_0
fparam = 1.454.*1e-4;
Rearth = 6378.1E3; % meters
alpha0=3.1E-4; %EOS - be careful because alpha is associated with frontogenesis too
beta0=7.4E-4;

xidxme=1:1420; % full zonal mean or subset of longitudes
nxm=length(xidxme)

RAC=RAC(xidxme,:);
hFacC=hFacC(xidxme,:,:);
RF=RF(xidxme,:,:);
DRF=DRF(xidxme,:,:);
RAS=RAS(xidxme,:);
RAW=RAW(xidxme,:);
hFacW=hFacW(xidxme,:,:);
hFacS=hFacS(xidxme,:,:);
DXC=DXC(xidxme,:);
DXG=DXG(xidxme,:);
DYC=DYC(xidxme,:);
DYG=DYG(xidxme,:);


% load mit coordinates
zmit=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','depth'));
ymit=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','latitude')).*111300;
latmit=repmat(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','latitude')',[nxm 1]);

lonmit=repmat(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','longitude'),[1 400]);
lonmit=lonmit(xidxme,:);

xmit=repmat(squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','longitude')).*111300, [1 400]);
xmit=xmit(xidxme,:);
xmit=xmit-mean(xmit(:));
xmit=xmit.*cosd(latmit);



%% load 3-D mitgcm U,V,W and generate mean budget tendencies <U><dU/dx>, <U><dV/dx>, <V><dU/dy> and <V><dV/dy>
% with < > representing a zonal and time mean

xmit3d=repmat(xmit,[1 1 136]);
ymit3d=repmat(ymit',[nxm 1  136]);

umit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','u'));
umit3d=umit3d(xidxme,:,:);
vmit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','v'));
vmit3d=vmit3d(xidxme,:,:);
wmit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','w'));
wmit3d=wmit3d(xidxme,:,:);

dvdy3d=zeros(size(vmit3d,1),size(vmit3d,2),size(vmit3d,3)+1);
dvdy3d(:,2:end-1,2:end)=(vmit3d(:,3:end,:)-vmit3d(:,1:end-2,:))./(ymit3d(:,3:end,:)-ymit3d(:,1:end-2,:));
dvdy3d(:,:,1)=dvdy3d(:,:,2);


dudxme=zeros(size(umit3d,1),size(umit3d,2),size(umit3d,3)+1);
dudxme(2:end-1,:,2:end)=(umit3d(3:end,:,:)-umit3d(1:end-2,:,:))./(xmit3d(3:end,:,:)-xmit3d(1:end-2,:,:));
zmitdwdz=cat(1,0,zmit);
dudxme(:,:,1)=dudxme(:,:,2);
wmex=cumtrapz(zmitdwdz,-dudxme,3);

 ududxme=squeeze(mean(umit3d,1,'omitnan')).*squeeze(mean(dudxme(:,:,2:end),1,'omitnan')); % < 3e-8
 dvdxme=zeros(size(umit3d));
 dvdxme(2:end-1,:,:)=(vmit3d(3:end,:,:)-vmit3d(1:end-2,:,:))./(xmit3d(3:end,:,:)-xmit3d(1:end-2,:,:));
 udvdxme=squeeze(mean(umit3d,'omitnan')).*squeeze(mean(dvdxme,1,'omitnan')); % <6e-9
% clear xmit3d dudxme dvdxme

dudyme=zeros(size(umit3d));
dudyme(:,2:end-1,:)=(umit3d(:,3:end,:)-umit3d(:,1:end-2,:))./(ymit3d(:,3:end,:)-ymit3d(:,1:end-2,:));
vdudyme=squeeze(mean(vmit3d,1,'omitnan')).*squeeze(mean(dudyme,1,'omitnan')); % < 2e-7

dvdyme=zeros(size(umit3d,1),size(umit3d,2),size(umit3d,3)+1);
dvdyme(:,2:end-1,2:end)=(vmit3d(:,3:end,:)-vmit3d(:,1:end-2,:))./(ymit3d(:,3:end,:)-ymit3d(:,1:end-2,:));
dvdyme(:,:,1)=dvdyme(:,:,2);
wmey=cumtrapz(zmitdwdz,-dvdyme,3);
vdvdyme=squeeze(mean(vmit3d,1,'omitnan')).*squeeze(mean(dvdyme(:,:,2:end),1,'omitnan')); % < 3e-8

% generate <W><dU/dz> and <W><dV/dz> mean budget terms

zmit3d=repmat(reshape(zmit,[1 1 136]),[nxm 400 1]);
dudzme=zeros(size(umit3d));
dudzme(:,:,2:end-1)=(umit3d(:,:,3:end)-umit3d(:,:,1:end-2))./(zmit3d(:,:,3:end)-zmit3d(:,:,1:end-2));
dudzme(:,:,1)=dudzme(:,:,2);
wdudzme=squeeze(mean(wmit3d,1,'omitnan')).*squeeze(mean(dudzme,1,'omitnan')); 
dvdzme=zeros(size(umit3d));
dvdzme(:,:,2:end-1)=(vmit3d(:,:,3:end)-vmit3d(:,:,1:end-2))./(zmit3d(:,:,3:end)-zmit3d(:,:,1:end-2));
wdvdzme=squeeze(mean(wmit3d,1,'omitnan')).*squeeze(mean(dvdzme,1,'omitnan')); % < 2e-8


%% generate zonal mean variables <U>, <V>, <T>, <W>, <S>
umit=squeeze(mean(umit3d,1,'omitnan'));
Tmit3d=ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','theta');
Tmit3d=Tmit3d(xidxme,:,:);
Tmit=squeeze(mean(Tmit3d,1,'omitnan'));

%% figure of MITgcm dV/dy for Appendix
figure('position',[10 10 600 350])
  contourf(latmit(1,:),zmit,squeeze(1e8.*mean(dvdyme(:,:,2:end),1,'omitnan'))',linspace(-50,50,21)); 
hold on
caxis([-50 50])
caxis(caxis)
hold on
[c,h]=contour(latmit(1,:),zmit,Tmit',-2:2:30,'linewidth',1,'color','k');
clabel(c,h,'fontsize',13);
hold on;
contour(latmit(1,:),zmit,umit',.2:.2:2,'linewidth',2,'color','w');
ylabel('depth (m)','FontSize',16,'FontName','Arial','FontWeight','normal');
xlabel('latitude ({\circ}N)','FontSize',16,'FontName','Arial','FontWeight','normal');
cbh = colorbar();
ylabel(cbh,'dv/dy (10^{-8} s^{-1})')
set(gca,'FontSize',16,'FontName','Arial','FontWeight','normal');
title('Zonal and time mean meridional divergence in the gcm','FontSize',16,'FontWeight','normal','FontName','Arial');
ylim([-300 0])
xlim([-8 8])
grid on
set(gcf,'color','w')
set(gcf,'PaperSize',[8 6]); 
set(gca,'colormap',cmocean('balance',20))
print(gcf,'MITgcm_dVdy.png','-dpng','-r200')



%% mitgcm and Eliassen model grids
y_gm=repmat(ymit',[length(zmit) 1]);
lat_gm=repmat(latmit(1,:),[length(zmit) 1]);
z_gm=repmat(zmit,[1 length(ymit)]);
% Eliassen grid is initially the same uniform horizontal grid with slight
% refinement in the vertical and then downsampled for faster matrix algebra
[y_g,z_g]=meshgrid(double(y_gm(1,:)),double([-0.75:-1.5:-600])');
[lat_g,z_g]=meshgrid(double(lat_gm(1,:)),double([-0.75:-1.5:-600])');

%%  zonal momentum diagnostics
% f and Coriolis force variables for
f_gm=fparam.*sind(lat_gm); % mitgcm
f=fparam.*sind(lat_g); % eliassen
vmit=squeeze(mean(vmit3d,1,'omitnan'));


% advection
Um_Advecwcor3d=squeeze(double(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','Um_Advec')));
Um_Advecwcor3d=Um_Advecwcor3d(xidxme,:,:);
Um_Cor3d=repmat(double(reshape(f_gm',[1 400 136])),[nxm 1 1]).*vmit3d;
%Um_Advecan3d=Um_Advecwcor-Um_Cor3d +repmat(reshape(vdudyme+wdudzme+ududxme,[1 400 136]),[nxm 1 1]);
Um_Cor=squeeze(mean(Um_Cor3d,1,'omitnan'));
Um_Advec=squeeze(mean(Um_Advecwcor3d,1,'omitnan'))-Um_Cor;
Um_Advecan=Um_Advec+vdudyme+wdudzme+ududxme;
TOTUTEND3d=squeeze(double(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','TOTUTEND'))./86400);
TOTUTEND3d=TOTUTEND3d(xidxme,:,:);

% vmix terms
VISRI_Um=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','VISrI_Um')));
VISRI_Um=VISRI_Um(xidxme,:,:);
VISRI_Um(:,:,1)=0;
Um_Impl=zeros(size(VISRI_Um));
Um_Impl(:,:,1:end-1)=(VISRI_Um(:,:,2:end)-VISRI_Um(:,:,1:end-1));
Um_Impl(:,:,end)=Um_Impl(:,:,end-1);% not bottom
normfac=double(squeeze(hFacW.*repmat(RAW,[1 1 136])).*DRF);
Um_Impl=Um_Impl./normfac;
Um_Ext=squeeze(double(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','Um_Ext')));
Um_Ext=Um_Ext(xidxme,:,:);
Um_Ext(isnan(Um_Ext))=0;

% dissipation terms (very small)
AB_gU=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','AB_gU')));
AB_gU=AB_gU(xidxme,:,:);
Um_Diss=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','Um_Diss')));
Um_Diss=Um_Diss(xidxme,:,:);

X=(squeeze(mean(Um_Impl+Um_Ext+Um_Diss+AB_gU,1,'omitnan'))+Um_Advecan)';
X(1,:)=X(2,:);
X0=squeeze(mean(Um_Impl+Um_Ext+Um_Diss+AB_gU,1,'omitnan'))';
X0(1,:)=X0(2,:);



% zonal PGF for vg
Um_dPHdx3d=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_ub.nc','Um_dPHdx')));
Um_dPHdx3d=Um_dPHdx3d(xidxme,:,:);
%ETAN=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_etan.nc','ETAN')));
%ETAN=ETAN(xidxme,:);
%Um_dPSdx3d=zeros(size(Um_dPHdx3d));
%Um_dPSdx3d(2:end-1,:,:)=-repmat(g.*(ETAN(3:end,:)-ETAN(1:end-2,:))./(2.*DXC(2:end-1,:)),[1 1 136]);
Um_dPsdx3dalt=-repmat(mean(squeeze(Um_Advecwcor3d+Um_dPHdx3d+...
    AB_gU+Um_Diss+Um_Ext+Um_Impl-TOTUTEND3d),3,'omitnan'),[1 1 136]);
Um_dPdx=squeeze(mean(Um_dPHdx3d+Um_dPsdx3dalt,1,'omitnan'))';
clear AB_gU Um_Diss Um_Ext Um_Impl Um_dPHdx3d Um_dPsdx3dalt TOTUTEND3d

%% meridional momentum budget terms

Vm_dPHdy3d=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','Vm_dPHdy')));
Vm_dPHdy3d=Vm_dPHdy3d(xidxme,:,:);

Vm_Diss=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','Vm_Diss')));
Vm_Diss=Vm_Diss(xidxme,:,:);


TOTVTEND3d=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','TOTVTEND')))./86400;
TOTVTEND3d=TOTVTEND3d(xidxme,:,:);

Vm_Advecwcor3d=squeeze(double(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','Vm_Advec')));
Vm_Advecwcor3d=Vm_Advecwcor3d(xidxme,:,:);
Vm_Cor3d=-repmat(double(reshape(f_gm',[1 400 136])),[nxm 1 1]).*umit3d;
Vm_Cor=squeeze(mean(Vm_Cor3d,1,'omitnan'));
Vm_Advec=squeeze(mean(Vm_Advecwcor3d,1,'omitnan'))-Vm_Cor;
Vm_Advecan=Vm_Advec+vdvdyme+wdvdzme+udvdxme;

VISRI_Vm=double(squeeze(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','VISrI_Vm')));
VISRI_Vm=VISRI_Vm(xidxme,:,:);
VISRI_Vm(:,:,1)=0;
Vm_Impl=zeros(size(VISRI_Vm));
Vm_Impl(:,:,1:end-1)=(VISRI_Vm(:,:,2:end)-VISRI_Vm(:,:,1:end-1));
Vm_Impl(:,:,end)=Vm_Impl(:,:,end-1);% not bottom
normfac=double(squeeze(hFacS.*repmat(RAS,[1 1 136])).*DRF);
Vm_Impl=Vm_Impl./normfac;
Vm_Ext=squeeze(double(ncread('~/OneDrive - NASA/documents/File_y1999_y2018_vb.nc','Vm_Ext')));
Vm_Ext=Vm_Ext(xidxme,:,:);
Vm_Ext(isnan(Vm_Ext))=0;



%AB_gV is inf for some reason but negligibly small
%Vm_dPSdy3d=zeros(size(Vm_dPHdy3d));
%Vm_dPSdy3d(:,2:end-1,:)=-repmat(g.*(ETAN(:,3:end)-ETAN(:,1:end-2))./(2.*DYC(:,2:end-1)),[1 1 136]);
Vm_dPsdy3dalt=-repmat(mean(squeeze(Vm_Advecwcor3d+Vm_dPHdy3d+...
    Vm_Diss+Vm_Ext+Vm_Impl-TOTVTEND3d),3,'omitnan'),[1 1 136]);
% there are very small mismatches between direct calculation and alt (residual
% calculation)... maybe issues with calculation of gradient from ETA?
% residual is chosen as more precise
Vm_dPdy=squeeze(mean(Vm_dPHdy3d+Vm_dPsdy3dalt,1,'omitnan'))';
Vm_Cor=Vm_Cor';

Y0=squeeze(mean(Vm_Impl+Vm_Ext+Vm_Diss,1,'omitnan'))';

Y=Y0+Vm_Advecan';

clear  Vm_Diss Vm_Ext Vm_Impl Vm_dPHdy3d Vm_dPSdy3d Vm_dPSdy3dalt TOTVTEND3d

%% thermal wind balance

TWBCor=zeros(size(Vm_Cor));
TWBCor(2:end-1,:)=(Vm_Cor(3:end,:)-Vm_Cor(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
TWBCor(1,:)=TWBCor(2,:);
TWBCor(end,:)=TWBCor(end,1,:);

TWBPGF=zeros(size(Vm_dPdy));
TWBPGF(2:end-1,:)=(Vm_dPdy(3:end,:)-Vm_dPdy(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
TWBPGF(1,:)=TWBPGF(2,:);
TWBPGF(end,:)=TWBPGF(end-1,:);

figure;
subplot(2,1,1),...
contourf(ymit./111300,zmit,-TWBCor,linspace(-2e-7,2e-7,81)); shading flat
hold on
caxis([-7e-8 7e-8])
caxis(caxis)
hold on
contour(ymit./111300,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
[c,h]=contour(ymit./111300,zmit,umit',linspace(-2,2,21),'linewidth',2,'color','k');
clabel(c,h);
ylabel('Depth (m)','FontSize',12);
xlabel('Latitude (^oN)','FontSize',12);
cbh = colorbar();
set(gca,'FontSize',12);
title('(a) time and zonal mean zonal f\partial u/\partial z (s^{-2})','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
colormap(cmocean('balance',28))


subplot(2,1,2),...
contourf(ymit./111300,zmit,-TWBPGF,linspace(-2e-7,2e-7,81)); shading flat
hold on
caxis([-7e-8 7e-8])
caxis(caxis)
hold on
contour(ymit./111300,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
[c,h]=contour(ymit./111300,zmit,umit',linspace(-2,2,21),'linewidth',2,'color','k');
clabel(c,h);
ylabel('Depth (m)','FontSize',12);
xlabel('Latitude (^oN)','FontSize',12);
cbh = colorbar();
set(gca,'FontSize',12);
title('(b) -\partial b/\partial y (s^{-2}) ','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
colormap(cmocean('balance',28))
set(gcf,'color','w')
print(gcf,'MITgcm_TWB.png','-dpng','-r200')

%% calculation of vg

vgmit=-(Um_dPdx)./f_gm;
vgmitf=zeros(size(vgmit));
vpoly=zeros(size(vgmit));
for i = 1:136
pvel=polyfit(latmit(700,[1:120,281:400]),vgmit(i,[1:120,281:400]),5);
vpoly(i,:)=polyval(pvel,latmit(700,:))';
vgmitf(i,[1:50,351:400])=vgmit(i,[1:50,351:400]);
wtp=((abs(ymit([51:140,261:350]))-abs(ymit(51)))./(abs(ymit(140))-abs(ymit(51))))';
wtg=1-wtp;
vgmitf(i,[51:140,261:350])=wtg.*squeeze(vgmit(i,[51:140,261:350]))+wtp.*squeeze(vpoly(i,[51:140,261:350]));
vgmitf(i,141:260)=vpoly(i,[141:260]);
vgmitf(i,:)=smooth(vgmitf(i,:),21);
end
dvgdy=zeros(size(vmit'));
dvgdy(:,2:end-1)=squeeze(vgmitf(:,3:end)-vgmitf(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
for i = 1:136
    dvgdy(i,:)=smooth(dvgdy(i,:)',11);
end
dvgdydwdz=zeros(137,400);
dvgdydwdz(2:end,:)=dvgdy;
dvgdydwdz(1,:)=dvgdy(1,:);
wgm=cumtrapz(zmitdwdz,-dvgdydwdz,1);
wgm=wgm(2:end,:);

vgdudyme=squeeze(vgmitf').*squeeze(mean(dudyme,1,'omitnan')); % < 2e-7

Xgm=-vgdudyme'-ududxme';


wmit=squeeze(mean(wmit3d,1,'omitnan'));


figure('position',[50 50 1000 500])
subplot(2,2,1),...
contourf(ymit./1000./111,zmit,vgmit,[-10,linspace(-0.1,0.1,21),10]); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
hold on
[c,h]=contour(ymit./1000./111,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
clabel(c,h);
contour(ymit./1000./111,zmit(2:end),86400.*squeeze(wmit(:,2:end))',linspace(-2,2,21),'linewidth',2,'color','k');
%clabel(c,h);
cbh = colorbar();

title('(a) meridional geostrophic velocity v_g (m/s)','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
colormap(cmocean('balance',40))
ylabel('depth (m)','FontSize',12);
xlabel('latitude (^oN)','FontSize',12);
set(gca,'FontSize',12);

subplot(2,2,2),...
contourf(ymit./1000./111,zmit,vgmitf,linspace(-0.1,0.1,21)); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
hold on
[c,h]=contour(ymit./1000./111,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
clabel(c,h);
contour(ymit./1000./111,zmit(2:end),86400.*squeeze(wmit(:,2:end))',linspace(-2,2,21),'linewidth',2,'color','k');
%clabel(c,h);
cbh = colorbar();
title('(b) v_g interpolated across the equator v_g^i (m/s)','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
colormap(cmocean('balance',40))
ylabel('depth (m)','FontSize',12);
xlabel('latitude (^oN)','FontSize',12);
set(gca,'FontSize',12);


subplot(2,2,4),...
contourf(ymit./1000./111,zmit,vmit'-vgmitf,linspace(-0.1,0.1,21)); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
hold on
[c,h]=contour(ymit./1000./111,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
clabel(c,h);
contour(ymit./1000./111,zmit(2:end),86400.*squeeze(wmit(:,2:end))',linspace(-2,2,21),'linewidth',2,'color','k');
%clabel(c,h);
cbh = colorbar();
title('(d) ageostrophic meridional velocity v-v_g^i (m/s)','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
ylabel('depth (m)','FontSize',12);
xlabel('latitude (^oN)','FontSize',12);
set(gca,'FontSize',12);



subplot(2,2,3),...
contourf(ymit./1000./111,zmit,vmit',linspace(-0.1,0.1,21)); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
hold on
[c,h]=contour(ymit./1000./111,zmit,Tmit',-2:2:30,'linewidth',1,'color','w');
clabel(c,h);
contour(ymit./1000./111,zmit(2:end),86400.*squeeze(wmit(:,2:end))',linspace(-2,2,21),'linewidth',2,'color','k');
cbh = colorbar();
title('(c) meridional velocity v (m/s)','FontSize',12,'FontWeight','normal');
ylim([-300 0])
xlim([-8 8])
grid on
set(gcf,'color','w')
ylabel('depth (m)','FontSize',12);
xlabel('latitude (^oN)','FontSize',12);
set(gca,'FontSize',12);
print(gcf,'MITgcm_interpolate_vg.png','-r200','-dpng')


%% temperature budget terms

TFLUX=double(ncread('~/OneDrive - NASA/documents/mitgcm_suf_avg.nc','TFLUX'));
TFLUX=TFLUX(xidxme,:);
surfQsw=double(ncread('~/OneDrive - NASA/documents/mitgcm_suf_avg.nc','oceQsw'));
surfQsw=surfQsw(xidxme,:);

ADVr_TH=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','ADVr_TH'));
ADVr_TH=ADVr_TH(xidxme,:,:);
ADVr_TH(:,:,1)=0;

ADVx_TH=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','ADVx_TH'));
ADVx_TH=ADVx_TH(xidxme,:,:);

ADVy_TH=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','ADVy_TH'));
ADVy_TH=ADVy_TH(xidxme,:,:);

CellVol = double(repmat(RAC,[1 1 136]).*DRF.*hFacC);

ADV_tend3d=zeros(size(CellVol));
ADV_tend3d(1:end-1,1:end-1,1:end-1)= ...
(ADVx_TH(2:end,1:end-1,1:end-1)-ADVx_TH(1:end-1,1:end-1,1:end-1))./CellVol(1:end-1,1:end-1,1:end-1) +...
(ADVy_TH(1:end-1,2:end,1:end-1)-ADVy_TH(1:end-1,1:end-1,1:end-1))./CellVol(1:end-1,1:end-1,1:end-1) +...
(ADVr_TH(1:end-1,1:end-1,1:end-1)-ADVr_TH(1:end-1,1:end-1,2:end))./CellVol(1:end-1,1:end-1,1:end-1);
ADV_tend3d=-ADV_tend3d;
ADV_tend3d(end,:,:)=ADV_tend3d(end-1,:,:);
ADV_tend3d(:,end,:)=ADV_tend3d(:,end-1,:);
ADV_tend3d(1,:,:)=ADV_tend3d(2,:,:);
ADV_tend3d(:,1,:)=ADV_tend3d(:,2,:);

ADV_tend=squeeze(mean(ADV_tend3d,1,'omitnan'))';


dTdyme=zeros(size(umit3d));
dTdxme=zeros(size(umit3d));
dTdyme(:,2:end-1,:)=(Tmit3d(:,3:end,:)-Tmit3d(:,1:end-2,:))./(ymit3d(:,3:end,:)-ymit3d(:,1:end-2,:));
dTdxme(2:end-1,:,:)=(Tmit3d(3:end,:,:)-Tmit3d(1:end-2,:,:))./(xmit3d(3:end,:,:)-xmit3d(1:end-2,:,:));
dTdyme(:,1,:)=dTdyme(:,2,:);
dTdyme(:,end,:)=dTdyme(:,end-1,:);
dTdxme(1,:,:)=dTdxme(2,:,:);
dTdxme(end,:,:)=dTdxme(end-1,:,:);

vdTdyme=squeeze(mean(vmit3d,1,'omitnan')).*squeeze(mean(dTdyme,1,'omitnan')); 
udTdxme=squeeze(mean(umit3d,1,'omitnan')).*squeeze(mean(dTdxme,1,'omitnan')); 

dTdzme=zeros(size(umit3d));
dTdzme(:,:,2:end-1)=(Tmit3d(:,:,3:end)-Tmit3d(:,:,1:end-2))./(zmit3d(:,:,3:end)-zmit3d(:,:,1:end-2));
dTdzme(:,:,1)=dTdzme(:,:,2);

wdTdzme=squeeze(mean(wmit3d,1,'omitnan')).*squeeze(mean(dTdzme,1,'omitnan')); %

DFrI_TH=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','DFrI_TH'));
DFrI_TH=DFrI_TH(xidxme,:,:);

KPPg_TH=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','KPPg_TH'));
KPPg_TH=KPPg_TH(xidxme,:,:);

TOTTTEND3d=double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','TOTTTEND'))./86400;
TOTTTEND3d=TOTTTEND3d(xidxme,:,:);

WTHMASS=squeeze(double(ncread('~/OneDrive - NASA/documents/mitgcm_hb_avg.nc','WTHMASS')));
WTHMASS=WTHMASS(xidxme,:,1);
surf_cor_tend=zeros(size(CellVol));
surf_cor_tend(:,:,1)=2.85e-9-WTHMASS./(DRF(:,:,1)); % the 2.85e-9 is just based on residual with TOTTTEND

DIF_tend3d=zeros(size(CellVol));
DFrI_TH(:,:,1)=0;
DIF_tend3d(1:end-1,1:end-1,1:end-1)=(DFrI_TH(1:end-1,1:end-1,1:end-1)-DFrI_TH(1:end-1,1:end-1,2:end))./CellVol(1:end-1,1:end-1,1:end-1);
DIF_tend3d(1:end-1,1:end-1,1:end-1)=-DIF_tend3d(1:end-1,1:end-1,1:end-1)...
    -(KPPg_TH(1:end-1,1:end-1,1:end-1)-KPPg_TH(1:end-1,1:end-1,2:end))./CellVol(1:end-1,1:end-1,1:end-1);
DIF_tend=squeeze(mean(DIF_tend3d,1,'omitnan'))';
DIF_tend(:,end)=DIF_tend(:,end-1);
DIF_tend3d(end,:,:)=DIF_tend3d(end-1,:,:);
DIF_tend3d(:,end,:)=DIF_tend3d(:,end-1,:);
%clear KPPg_TH DFrI_TH WTHMASS

swfrac = 0.62.*exp(double(RF(:,:,1:end-1))./0.6) + (1.0-0.62).*exp(double(RF(:,:,1:end-1))./20.0);
swfrac1 = 0.62.*exp(double(RF(:,:,2:end))./0.6) + (1.0-0.62).*exp(double(RF(:,:,2:end))./20.0);


Qsw_tend3d=zeros(size(CellVol));
Qsw_tend3d(:,:,1:end-1)=repmat(surfQsw,[1 1 135])./3994./rhoref./DRF(:,:,1:end-1)./hFacC(:,:,1:end-1).*(swfrac-swfrac1);
Tflx_tend3d=zeros(size(CellVol));
Tflx_tend3d(:,:,1)=(TFLUX-surfQsw)./rhoref./3994./DRF(:,:,1)./hFacC(:,:,1);


%clear ADV_tendan3d
%clear RHS3d swfrac swfrac1 hFacC DRF

% salt offline
load('/Users/dbwhitt/OneDrive - NASA/documents/mitgcm_offline_salt_flux_avg_anoms.mat');
SFnow=SFnow(xidxme,:,:,:);
Smit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','salt'));
Smit3d=Smit3d(xidxme,:,:);
Smit=squeeze(mean(Smit3d,1,'omitnan'));

dSdyme=zeros(size(umit3d));
dSdxme=zeros(size(umit3d));
dSdyme(:,2:end-1,:)=(Smit3d(:,3:end,:)-Smit3d(:,1:end-2,:))./(ymit3d(:,3:end,:)-ymit3d(:,1:end-2,:));
dSdxme(2:end-1,:,:)=(Smit3d(3:end,:,:)-Smit3d(1:end-2,:,:))./(xmit3d(3:end,:,:)-xmit3d(1:end-2,:,:));
dSdyme(:,1,:)=dSdyme(:,2,:);
dSdyme(:,end,:)=dSdyme(:,end-1,:);
dSdxme(1,:,:)=dSdxme(2,:,:);
dSdxme(end,:,:)=dSdxme(end-1,:,:);


vdSdyme=squeeze(mean(vmit3d,1,'omitnan')).*squeeze(mean(dSdyme,1,'omitnan')); 
udSdxme=squeeze(mean(umit3d,1,'omitnan')).*squeeze(mean(dSdxme,1,'omitnan')); 


dSdzme=zeros(size(umit3d));
dSdzme(:,:,2:end-1)=(Smit3d(:,:,3:end)-Smit3d(:,:,1:end-2))./(zmit3d(:,:,3:end)-zmit3d(:,:,1:end-2));
dSdzme(:,:,1)=dSdzme(:,:,2);
wdSdzme=squeeze(mean(wmit3d,1,'omitnan')).*squeeze(mean(dSdzme,1,'omitnan')); %
ADV_tend_salt3d=zeros(size(ADV_tend3d));
ADV_tend_salt3d(2:end-1,2:end-1,2:end-1)=-(...
    (SFnow(3:end,2:end-1,2:end-1,1)-SFnow(1:end-2,2:end-1,2:end-1,1))./(xmit3d(3:end,2:end-1,2:end-1)-xmit3d(1:end-2,2:end-1,2:end-1)) +...
    (SFnow(2:end-1,3:end,2:end-1,2)-SFnow(2:end-1,1:end-2,2:end-1,2))./(ymit3d(2:end-1,3:end,2:end-1)-ymit3d(2:end-1,1:end-2,2:end-1)) + ...
    (SFnow(2:end-1,2:end-1,3:end,3)-SFnow(2:end-1,2:end-1,1:end-2,3))./(zmit3d(2:end-1,2:end-1,3:end)-zmit3d(2:end-1,2:end-1,1:end-2)) ); 
ADV_tend_salt3d(:,:,1)=ADV_tend_salt3d(:,:,2);
ADV_tend_salt2d=squeeze(mean(ADV_tend_salt3d,1,'omitnan'));
clear umit3d wmit3d vmit3d SFnow 
DIA_salt3d=-ADV_tend_salt3d;

for i = 1:nxm
galpha3d(i,:,:)=g.*sw_alpha(squeeze(Smit3d(i,:,:)),squeeze(Tmit3d(i,:,:)),zeros(400,136));
gbeta3d(i,:,:)=g.*sw_beta(squeeze(Smit3d(i,:,:)),squeeze(Tmit3d(i,:,:)),zeros(400,136));
end
galpha2d=squeeze(mean(galpha3d,1,'omitnan'));
gbeta2d=squeeze(mean(gbeta3d,1,'omitnan'));

Bb0=(galpha2d').*squeeze(mean((Tflx_tend3d+Qsw_tend3d+DIF_tend3d+surf_cor_tend),1,'omitnan'))'-(gbeta2d').*squeeze(mean(DIA_salt3d,1,'omitnan'))';

clear Tflx_tend3d Qsw_tend3d surf_cor_tend DIF_tend3d

Bb=Bb0+(galpha2d').*squeeze(mean(ADV_tend3d,1,'omitnan'))'+(galpha2d.*wdTdzme)'+(galpha2d.*udTdxme)'+(galpha2d.*vdTdyme)'-...
        (gbeta2d').*squeeze(mean(ADV_tend_salt3d,1,'omitnan'))'-(gbeta2d.*udSdxme)'-(gbeta2d.*vdSdyme)'-(gbeta2d.*wdSdzme)';

% mostly Bgmx:
Bgm=+(gbeta2d.*udSdxme)'...
    -(galpha2d.*udTdxme)'...
     +(gbeta2d').*squeeze(vgmitf).*(squeeze(mean(dSdyme,1,'omitnan'))')...
     -(galpha2d').*squeeze(vgmitf).*(squeeze(mean(dTdyme,1,'omitnan'))');

%% smooth all forcings at the tendency equation level with a 1/2 deg moving average in latitude
for i = 1:136
    Bb(i,:)=smooth(squeeze(double(Bb(i,:))),11);
    Bb(i,1:6)=Bb(i,7);
    Bb(i,end-6:end)=Bb(i,end-7);
    Bb0(i,:)=smooth(squeeze(double(Bb0(i,:))),11);
    Bb0(i,1:6)=Bb0(i,7);
    Bb0(i,end-6:end)=Bb0(i,end-7);
    Bgm(i,:)=smooth(squeeze(double(Bgm(i,:))),11);
    Bgm(i,1:6)=Bgm(i,7);
    Bgm(i,end-6:end)=Bgm(i,end-7);
    X(i,:)=smooth(squeeze(double(X(i,:))),11);
    X0(i,:)=smooth(squeeze(double(X0(i,:))),11);
    Xgm(i,:)=smooth(squeeze(double(Xgm(i,:))),11);
    Y(i,:)=smooth(squeeze(double(Y(i,:))),11);
    Y0(i,:)=smooth(squeeze(double(Y0(i,:))),11);
end

%% define and smooth Eliassen operator coefficients S2,F2,N2
buoy3d=zeros(size(Tmit3d));
for k=1:136
buoy3d(:,:,k)=-9.81.*sw_dens0(squeeze(Smit3d(:,:,k)),squeeze(Tmit3d(:,:,k)))./rhoref;
end
buoy=squeeze(mean(-9.81.*sw_dens0(Smit3d,Tmit3d)./rhoref,1,'omitnan'));
buoy=buoy';
umit=umit';
for i = 1:136
    umit(i,:)=smooth(squeeze(double(umit(i,:))),11);
    buoy(i,:)=smooth(squeeze(double(buoy(i,:))),11);
end
dbdy=zeros(size(buoy));
dudz=zeros(size(buoy));
dudy=zeros(size(buoy));
dbdz=zeros(size(buoy));

dbdy(:,2:end-1)=(buoy(:,3:end)-buoy(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
dudz(2:end-1,:)=(umit(3:end,:)-umit(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
dudy(:,2:end-1)=(umit(:,3:end)-umit(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
dbdz(2:end-1,:)=(buoy(3:end,:)-buoy(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
dbdz(1,:)=dbdz(2,:);
dbdz(end,:)=dbdz(end-1,:);

DzY=zeros(size(Y));
DzY(2:end-1,:)=(Y(3:end,:)-Y(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
DzY(1,:)=DzY(2,:)+ DzY(2,:)-DzY(3,:);
DzY(end,:)=DzY(end-1,:)+ DzY(end-1,:)-DzY(end-2,:);

DzY0=zeros(size(Y));
DzY0(2:end-1,:)=(Y0(3:end,:)-Y0(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
DzY0(1,:)=DzY0(2,:)+ DzY0(2,:)-DzY0(3,:);
DzY0(end,:)=DzY0(end-1,:)+ DzY0(end-1,:)-DzY0(end-2,:);


% these are defined assuming turbulent thermal wind balance rather than
% thermal wind balance
S2fb=interp2(double(y_gm),double(z_gm),0.5.*double(DzY0-2.*dbdy),double(y_g),double(z_g));
S2fb(1,:)=S2fb(2,:);
S2fb(end,:)=S2fb(end-1,:);

S2fu=f.*interp2(double(y_gm),double(z_gm),double(dudz),double(y_g),double(z_g))-...
    0.5.*interp2(double(y_gm),double(z_gm),double(DzY0),double(y_g),double(z_g));
S2fu(1,:)=S2fu(2,:);
S2fu(end,:)=S2fu(end-1,:);

S2b=interp2(double(y_gm),double(z_gm),double(-dbdy),double(y_g),double(z_g));
S2b(1,:)=S2b(2,:);
S2u=f.*interp2(double(y_gm),double(z_gm),double(dudz),double(y_g),double(z_g));
S2u(1:2,:)=repmat(S2u(3,:),[2 1]);
F2=f.*(f-interp2(double(y_gm),double(z_gm),double(dudy),double(y_g),double(z_g)));
F2(1,:)=F2(2,:);
N2=interp2(double(y_gm),double(z_gm),double(dbdz),double(y_g),double(z_g));
testpdf=normpdf(y_g,0,3e5)./max(normpdf(y_g(:),0,3e5));
S2=S2u.*(testpdf)+S2b.*(1-testpdf);
S2(1,:)=S2(2,:);
S2(end,:)=S2(end-1,:);


dy=y_g(1,2)-y_g(1,1);
dz=z_g(2,1)-z_g(1,1);


%PV=F2.*N2-S2.^2;


%% define and interpolate forcing operators 


DzX=zeros(size(X));
DzX(2:end-1,:)=(X(3:end,:)-X(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
DzX(1,:)=DzX(2,:)+ DzX(2,:)-DzX(3,:);
DzX(end,:)=DzX(end-1,:)+ DzX(end-1,:)-DzX(end-2,:);

DzX0=zeros(size(X));
DzX0(2:end-1,:)=(X0(3:end,:)-X0(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
DzX0(1,:)=DzX0(2,:)+ DzX0(2,:)-DzX0(3,:);
DzX0(end,:)=DzX0(end-1,:)+ DzX0(end-1,:)-DzX0(end-2,:);

DzXgm=zeros(size(X));
DzXgm(2:end-1,:)=(Xgm(3:end,:)-Xgm(1:end-2,:))./(z_gm(3:end,:)-z_gm(1:end-2,:));
DzXgm(1,:)=DzXgm(2,:)+ DzXgm(2,:)-DzXgm(3,:);
DzXgm(end,:)=DzXgm(end-1,:)+ DzXgm(end-1,:)-DzXgm(end-2,:);

DyB=zeros(size(X));
DyB(:,2:end-1)=(Bb(:,3:end)-Bb(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
DyB(:,1)=DyB(:,2);
DyB(:,end)=DyB(:,end-1);

DyB0=zeros(size(X));
DyB0(:,2:end-1)=(Bb0(:,3:end)-Bb0(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
DyB0(:,1)=DyB0(:,2);
DyB0(:,end)=DyB0(:,end-1);

DyBgm=zeros(size(Bgm));
DyBgm(:,2:end-1)=(Bgm(:,3:end)-Bgm(:,1:end-2))./(y_gm(:,3:end)-y_gm(:,1:end-2));
DyBgm(:,1)=DyBgm(:,2);
DyBgm(:,end)=DyBgm(:,end-1);

DzX=interp2(double(y_gm),double(z_gm),double(DzX),double(y_g),double(z_g));
DzX(1:4,:)=repmat(DzX(4,:),[4 1]);
DzX0=interp2(double(y_gm),double(z_gm),double(DzX0),double(y_g),double(z_g));
DzX0(1:4,:)=repmat(DzX0(4,:),[4 1]);
DzXgm=interp2(double(y_gm),double(z_gm),double(DzXgm),double(y_g),double(z_g));
DzXgm(1:4,:)=repmat(DzXgm(4,:),[4 1]);

% DzY=interp2(double(y_gm),double(z_gm),double(DzY),double(y_g),double(z_g));
% DzY(1:4,:)=repmat(DzY(4,:),[4 1]);
% DzY0=interp2(double(y_gm),double(z_gm),double(DzY0),double(y_g),double(z_g));
% DzY0(1:4,:)=repmat(DzY0(4,:),[4 1]);


DyB=interp2(double(y_gm),double(z_gm),double(DyB),double(y_g),double(z_g));
DyB(1,:)=DyB(2,:);
DyB(:,end-1)=DyB(:,end-2);
DyB(:,end)=DyB(:,end-1);

DyBgm=interp2(double(y_gm),double(z_gm),double(DyBgm),double(y_g),double(z_g));
DyBgm(1,:)=DyBgm(2,:);
DyBgm(:,end-1)=DyBgm(:,end-2);
DyBgm(:,end)=DyBgm(:,end-1);

DyB0=interp2(double(y_gm),double(z_gm),double(DyB0),double(y_g),double(z_g));
DyB0(1,:)=DyB0(2,:);
DyB0(:,end-1)=DyB0(:,end-2);
DyB0(:,end)=DyB0(:,end-1);


%% set viscosity m2/s for the Eliassen model
nuz=1e-4.*ones(size(umit)); % must be constant for now
nuz=interp2(double(y_gm),double(z_gm),double(nuz),double(y_g),double(z_g));


%% script to 2olve Forced Eliassen-Sawyer Problem in arbitrary background
% flow as in Whitt and Thomas (2013) JPO
% Dan Whitt 

 
Iplt=2:2:400;
Jplt=2:2:400;

f=f(Jplt,Iplt);
F2 = F2(Jplt,Iplt);
S2 = S2(Jplt,Iplt);
N2 = N2(Jplt,Iplt);
nuz=nuz(Jplt,Iplt);
y_g = y_g(Jplt,Iplt);
z_g = z_g(Jplt,Iplt);
dz = z_g(2,1)-z_g(1,1)
dy = y_g(1,2)-y_g(1,1)
szY = size(F2,2)
szZ=size(F2,1)

DzX=DzX(Iplt,Jplt);
DyB =DyB(Iplt,Jplt);

DzX0=DzX0(Iplt,Jplt);
DyB0 =DyB0(Iplt,Jplt);

DzXgm=DzXgm(Iplt,Jplt);
DyBgm =DyBgm(Iplt,Jplt);




% operators
A = [-1.*ones(szY,1),16.*ones(szY,1),-30.*ones(szY,1),16.*ones(szY,1),-1.*ones(szY,1)];
K2y = (1./(12.*dy.^2)).*spdiags(A,[-2,-1,0,1,2],szY,szY);

% dpsidy=-w=0 at the left and right side boundaries -> Neumann condition
K2y(1, 1) = K2y(1,1)+K2y(1,2);
K2y(1, 2) = K2y(1,2)+K2y(1,3);
K2y(2,1) =  K2y(2,1)+ K2y(2,4);

K2y(end, end) = K2y(end,end)+K2y(end,end-1);
K2y(end, end-1) = K2y(end,end-1)+K2y(end,end-2);
K2y(end-1, end) = K2y(end-1,end)+K2y(end-1,end-3); 


B = [ones(szY,1),-8.*ones(szY,1),zeros(szY,1),8.*ones(szY,1),-1.*ones(szY,1)];
K1y = (1./(12.*dy)).*spdiags(B,[-2,-1,0,1,2],szY,szY);

% dpsidy=-w=0 at the left and right side boundaries -> Neumann condition
K1y(1, 1) = K1y(1,1)-K1y(1,2);
K1y(1, 2) = K1y(1,2)-K1y(1,3);
K1y(2,1) =  K1y(2,1)- K1y(2,4);

K1y(end, end) = K1y(end,end)-K1y(end,end-1);
K1y(end, end-1) = K1y(end,end-1)-K1y(end,end-2);
K1y(end-1, end) = K1y(end-1,end)-K1y(end-1,end-3); 

Dy = kron(speye(szZ),K1y);
Dyy = kron(speye(szZ),K2y);

B = fliplr(B);

K1z = (1./(12.*dz)).*spdiags(B,[-2,-1,0,1,2],szZ,szZ);
K2z = (1./(12.*dz.^2)).*spdiags(A,[-2,-1,0,1,2],szZ,szZ);
clear A
% Dirichlet: psi = 0 = constant (dpsi/dy =-w=0) at the top boundary
% no changes needed

% Neumann:  dpsi/dz = v = 0 at the bottom boundary

K2z(1, 1) = K2z(1,1)+K2z(1,2);
K2z(1, 2) = K2z(1,2)+K2z(1,3);
K2z(2,1) =  K2z(2,1)+ K2z(2,4);
K1z(1, 1) = K1z(1,1)-K1z(1,2);
K1z(1, 2) = K1z(1,2)-K1z(1,3);
K1z(2,1) =  K1z(2,1)- K1z(2,4);


Dzz = kron(K2z,speye(szY));
Dz = kron(K1z,speye(szY));
N2vec = vecES(N2);
S2vec = vecES(S2);
F2vec = vecES(F2);
yvec=vecES(y_g);

%RHS of E-S equation
b = vecES(-f.*(DzX) - (DyB));
bE = vecES(-f.*(DzX-DzX0) - (DyB-DyB0));
bEk = vecES(-f.*DzX0);
bEddyMom=vecES(-f.*(DzX-DzX0));
bEB=vecES(-DyB+DyB0);
bB0=vecES(-DyB0);
DzX0m=repmat(mean(DzX0,2),[1 200]);
bEknoy = vecES(-f.*DzX0m);
bg=vecES(-f.*(DzXgm)-(DyBgm));

% Elliassen-Sawyer Operator
A=sparse(bsxfun(@times,F2vec,Dzz) + 2.*bsxfun(@times,S2vec,Dz*Dy) + bsxfun(@times,N2vec,Dyy));

% Add beta plane term
% approximation beta ~ f/y
A=A+sparse(bsxfun(@times,S2vec./yvec,Dz));

% add viscosity to damp out small scales in secondary circulation (viscosity only operates on perturbation hydrostatic zonal vorticity omega_x=dv/dz, it does not operate on perturbation zonal momentum or buoyancy) 
nuzvec = vecES(nuz);
% add vertical frictional operator to equilibrate the secondary circulation scales
A = A + sparse(bsxfun(@times,nuzvec.^2,Dzz*Dzz*Dzz));



display('Solve Ax=b exactly')
psivec = A\b;
psivecEk=A\bEk;
psivecEknoy=A\bEknoy;
psivecEM=A\bEddyMom;
psivecB0=A\bB0;
psivecEB=A\bEB;
psivecE=A\bE;
psivecg = A\bg;


clear A
psi = matESH(psivec);
dpsidzvec = Dz*psivec;
dpsidyvec = Dy*psivec;
w = -dpsidyvec; % Whitt and Thomas 2013 signconvention 
v = dpsidzvec;
wEk=-Dy*psivecEk;
wEknoy=-Dy*psivecEknoy;
dvdyEk=Dy*(Dz*psivecEk);
dvdyEknoy=Dy*(Dz*psivecEknoy);


wEM=-Dy*psivecEM;
wB0=-Dy*psivecB0;
wEB=-Dy*psivecEB;
wE=-Dy*psivecE;
wg=-Dy*psivecg;

vg=Dz*psivecg;

wmexg=wmex;
wmexg(abs(86400.*wmexg)>5)=nan; % the threshold eliminates some noise in the zonal means resulting from strong vertical motion at islands
wmexg=squeeze(mean(wmexg,1,'omitnan'))';
wmeyg=wmey; % the threshold eliminates some noise in the zonal means resulting from strong vertical motion at islands
wmeyg(abs(86400.*wmeyg)>5)=nan;
wmeyg=squeeze(mean(wmeyg,1,'omitnan'))';


%pause
%% plot
figure('position',[50 50 1000 500])
subplot(2,2,1),...
contourf(y_g./111300,z_g,real(matESH(v)),linspace(-.3,.3,61),'linestyle','none'); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
hold on
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
ylim([-300 0])
xlim([-8 8])
grid on
title('(a) Eliassen v_a (m/s)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))
hold on;
caxis(caxis);
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(linspace(-8,8,4),-50.*ones(4,1),'k-','linewidth',2)


subplot(2,2,2),...
contourf(y_g./111300,z_g,real(matESH(w).*86400),linspace(-3,3,31),'linestyle','none'); shading flat
hold on
caxis([-2 2])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('(b) Eliassen w_a (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))
hold on;
caxis(caxis);
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(linspace(-8,8,4),-50.*ones(4,1),'k-','linewidth',2)


subplot(2,2,3),...
contourf(y_gm./111300,z_gm,vmit'-vgmitf,linspace(-.3,.3,61),'linestyle','none'); shading flat
hold on
caxis([-.1 .1])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
ylim([-300 0])
xlim([-8 8])
grid on
title('(c) MITgcm v_a (m/s)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))
hold on;
caxis(caxis);
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(linspace(-8,8,4),-50.*ones(4,1),'k-','linewidth',2)

subplot(2,2,4),...
contourf(y_gm./111300,z_gm,((wmeyg(2:end,:)-wgm).*86400),linspace(-3,3,31),'linestyle','none'); shading flat
hold on
caxis([-2 2])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');

cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
hold on;
caxis(caxis);
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(linspace(-8,8,4),-50.*ones(4,1),'k-','linewidth',2)

colormap(gca,cmocean('balance',20))
title('(d) MITgcm w_a m/d','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
set(gcf,'color','w')
print(gcf,'EliassenvsMITgcm_vawa.png','-r200','-dpng')

figure('position',[50 50 600 300]);
wESout=real(matESH(w).*86400);
plot(y_g(10,101:200)./111500,wESout(17,101:200),'b-','LineWidth',4);
 hold on
 plot(y_gm(10,201:400)./111500,86400.*(squeeze(wmeyg(20,201:400)-wgm(20,201:400))),'b-','LineWidth',1);
 plot(-y_g(10,1:100)./111500,wESout(17,1:100),'r--','LineWidth',4);
 hold on
 plot(-y_gm(10,1:200)./111500,86400.*(squeeze(wmeyg(20,1:200)-wgm(20,1:200))),'r--','LineWidth',1);
ylabel('w_a at 50 m (m/d)')
 xlabel('Latitude from the equator')
 legend('Northern hemisphere - Eliassen','Northern hemisphere - MITgcm','Southern hemisphere - Eliassen','Southern hemisphere - MITgcm')
set(gcf,'color','w')
xlim([0 8])
set(gca,'fontsize',14)
grid on
% Create doublearrow
annotation(gcf,'doublearrow',[0.321666666666667 0.32],...
    [0.733333333333334 0.45],'LineWidth',3,'Head2Width',15,'Head2Length',15,...
    'Head1Width',15,...
    'Head1Length',15);

% Create textarrow
annotation(gcf,'textarrow',[0.440000000000001 0.336666666666668],...
    [0.590000000000002 0.553333333333336],...
    'String',{'Cross-equatorial','asymmetry in upwelling'},...
    'FontSize',18);
print(gcf,'Eliassen_evaluation_lineplot.png','-r200','-dpng')

 figure('position',[10 10 600 900])
 subplot(2,1,1),...
 plot(y_g(10,101:200)./111500,wESout(17,101:200),'b-','LineWidth',4);
 hold on
 plot(-y_g(10,1:100)./111500,wESout(17,1:100),'r-','LineWidth',4);

wESEk=real(matESH(wEk).*86400);
plot(y_g(10,101:200)./111500,wESEk(17,101:200),'b:','LineWidth',2);
 hold on
 plot(-y_g(10,1:100)./111500,wESEk(17,1:100),'r:','LineWidth',2);

 wESE=real(matESH(wE).*86400);
 wESB0=real(matESH(wB0).*86400);

plot(y_g(10,101:200)./111500,wESE(17,101:200),'b-.','LineWidth',2);
 hold on
 plot(-y_g(10,1:100)./111500,wESE(17,1:100),'r-.','LineWidth',2);

 ylim([-1 2])
 ylabel('w_a at 50 m (m/d)')
 xlabel('Latitude from Equator')
 title('(a) Eliassen vertical velocity at 50 m in each hemisphere','fontweight','normal')
xlim([0 6])
grid on
set(gcf,'color','w')
set(gca,'Fontsize',15)

set(gca,'fontsize',14)

 subplot(2,1,2),...
     plot(y_g(10,101:200)./111500,wESout(17,101:200)-wESout(17,100:-1:1),'k-','LineWidth',4);
 hold on
plot(y_g(10,101:200)./111500,wESEk(17,101:200)-wESEk(17,100:-1:1),'k:','LineWidth',2);
 hold on
plot(y_g(10,101:200)./111500,wESE(17,101:200)-wESE(17,100:-1:1),'k-.','LineWidth',2);

 ylim([-1 2])
 ylabel('w_a at 50 m (m/d)')
 xlabel('Latitude from Equator')
 title('(b) Cross-equator asymmetry (north minus south)', 'fontweight','normal')

legend1=legend('All drivers','Wind: X_{vmix}','Eddies: X_{eddy} and B_{eddy}')

grid on
xlim([0 6])
%wESB0=real(matESH(wB0).*86400);
%plot(y_g(10,101:200)./111500,wESB0(17,101:200),'green-.','LineWidth',2);
% hold on
% plot(-y_g(10,1:100)./111500,wESB0(17,1:100),'yellow-.','LineWidth',2);
set(gca,'fontsize',14)

set(legend1,'Position',[0.587500000000002 0.7725 0.294166666666667 0.07625]);

% Create textarrow
annotation(gcf,'textarrow',[0.298333333333334 0.346666666666667],...
    [0.63625 0.705],'String',{'Eddy-driven','asymmetry'},'FontSize',18,...
    'FontName','Helvetica Neue');

% Create doublearrow
annotation(gcf,'doublearrow',[0.356666666666668 0.356666666666667],...
    [0.7775 0.645],'LineWidth',3,'Head2Width',15,'Head2Length',15,...
    'Head1Width',15,...
    'Head1Length',15);

% Create textbox
annotation(gcf,'textbox',...
    [0.572666666666672 0.84 0.152333333333333 0.0412500000000005],...
    'Color',[1 0 0],...
    'String','Red=South',...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.739333333333339 0.8325 0.152333333333333 0.0487500000000004],...
    'Color',[0 0 1],...
    'String','Blue=North',...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textarrow
annotation(gcf,'textarrow',[0.378333333333333 0.303333333333333],...
    [0.87125 0.82],'String',{'Wind/background','asymmetry'},'FontSize',18,...
    'FontName','Helvetica Neue');

% Create doublearrow
annotation(gcf,'doublearrow',[0.295000000000002 0.298333333333334],...
    [0.84191666666667 0.794583333333335],'LineWidth',3,'Head2Width',15,...
    'Head2Length',15,...
    'Head1Width',15,...
    'Head1Length',15);
set(gca,'Fontsize',15)
print(gcf,'Eliassen_thesis.png','-r200','-dpng')

%print(gcf,'Eliassen_eastvswest.png','-r200','-dpng')




figure('position',[40 40 1200 500]);
subplot(2,3,1),...
contourf(y_g./111300,z_g,real(matESH(wEM+wEB).*86400),linspace(-3,3,61),'linestyle','none'); shading flat
hold on
caxis([-1 1])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
grid on
contour(y_g./111300,z_g,real(matESH(psivecE)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecE)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
title('(a) Eliassen w_a [Eddy fluxes] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,4),...
contourf(y_g./111300,z_g,-f.*(DzX-DzX0)-(DyB-DyB0),linspace(-2e-13,2e-13,51),'linestyle','none'); shading flat
hold on
caxis([-8e-14 8e-14])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecE)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecE)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
grid on
title('(d) -f \partial_z X_{eddy} -\partial_y B_{eddy} (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,2),...
contourf(y_g./111300,z_g,real(matESH(wEB).*86400),linspace(-3,3,61),'linestyle','none'); shading flat
hold on
caxis([-1 1])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecEB)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEB)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
grid on
title('(c) Eliassen w_a [Eddy fluxes buoyancy] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,5),...
contourf(y_g./111300,z_g,-(DyB-DyB0),linspace(-2e-13,2e-13,51),'linestyle','none'); shading flat
hold on
caxis([-8e-14 8e-14])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
contour(y_g./111300,z_g,real(matESH(psivecEB)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEB)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('(e)  -\partial_y B_{eddy} (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))


subplot(2,3,3),...
contourf(y_g./111300,z_g,real(matESH(wEM).*86400),linspace(-3,3,61),'linestyle','none'); shading flat
hold on
caxis([-1 1])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
contour(y_g./111300,z_g,real(matESH(psivecEM)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEM)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('(c) Eliassen w_a [Eddy fluxes momentum] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,6),...
contourf(y_g./111300,z_g,-f.*(DzX-DzX0),linspace(-2e-13,2e-13,51),'linestyle','none'); shading flat
hold on
caxis([-8e-14 8e-14])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
%contour(y_g./111300,z_g,real(matESH(psivecEM)),linspace(-2,2,21),'linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEM)),linspace(-5,-0.2,25),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEM)),linspace(0,5,26),'-','linewidth',1,'color','w'); 
grid on
title('(f)  -f\partial_z X_{eddy} (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))
set(gcf,'color','w')
print(gcf,'Eliassen_eddyforcings.png','-r200','-dpng')








figure('position',[40 40 1300 500]);
subplot(2,3,1),...
contourf(y_g./111300,z_g,real(matESH(wEk).*86400),linspace(-3,3,31),'linestyle','none'); shading flat
hold on
caxis([-2 2])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecEk)),linspace(-5,-1,5),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEk)),linspace(0,5,6),'-','linewidth',1,'color','w'); 

grid on
title('(a) Eliassen w_a [Wind/Vertical Mixing Momentum] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,4),...
contourf(y_g./111300,z_g,-f.*(DzX0),linspace(-1e-12,1e-12,51),'linestyle','none'); shading flat
hold on
caxis([-4e-13 4e-13])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecEk)),linspace(-5,-1,5),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEk)),linspace(0,5,6),'-','linewidth',1,'color','w'); 

grid on
title('(d) -f \partial_z X_{vmix} (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))



subplot(2,3,2),...
contourf(y_g./111300,z_g,real(matESH(wB0).*86400),linspace(-3,3,61),'linestyle','none'); shading flat
hold on
caxis([-1 1])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.25:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecB0)),linspace(-5,-0.1,50),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecB0)),linspace(0,5,51),'-','linewidth',1,'color','w'); 

grid on
title('(b) Eliassen w_a [Vertical mixing buoyancy] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,5),...
contourf(y_g./111300,z_g,-(DyB0),linspace(-1e-13,1e-13,51),'linestyle','none'); shading flat
hold on
caxis([-4e-14 4e-14])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.25:29,'linewidth',1,'color','k');
cbh = colorbar();
contour(y_g./111300,z_g,real(matESH(psivecB0)),linspace(-5,-0.1,50),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecB0)),linspace(0,5,51),'-','linewidth',1,'color','w'); 

grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('(e)  -\partial_y B_{vmix} (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))


subplot(2,3,3),...
contourf(y_g./111300,z_g,real(matESH(wEknoy).*86400),linspace(-4,4,81),'linestyle','none'); shading flat
hold on
caxis([-3 3])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
contour(y_g./111300,z_g,real(matESH(psivecEknoy)),linspace(-5,-1,5),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEknoy)),linspace(0,5,6),'-','linewidth',1,'color','w'); 

grid on
title('(c) Eliassen w_a [horizontally-uniform X_{vmix}] (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))

subplot(2,3,6),...
contourf(y_g./111300,z_g,-f.*(DzX0m),linspace(-1e-12,1e-12,51),'linestyle','none'); shading flat
hold on
caxis([-4e-13 4e-13])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.5:29,'linewidth',1,'color','k');
cbh = colorbar();
contour(y_g./111300,z_g,real(matESH(psivecEknoy)),linspace(-5,-1,5),'--','linewidth',1,'color','w'); 
contour(y_g./111300,z_g,real(matESH(psivecEknoy)),linspace(0,5,6),'-','linewidth',1,'color','w'); 
grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('(f)  -f<\partial_z X_{vmix}>_y (s^{-3})','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))


set(gcf,'color','w')
print(gcf,'Eliassen_vmixforcings_6panel.png','-r200','-dpng')




figure;
contourf(y_gm./111300,z_gm,(wmit.*86400)',linspace(-3,3,31),'linestyle','none'); shading flat
hold on
caxis([-2 2])
caxis(caxis)
contour(y_gm./111300,z_gm,-rhoref.*buoy./g -1000,20:.25:29,'linewidth',1,'color','k');
cbh = colorbar();
grid on
ylim([-300 0])
xlim([-8 8])
grid on
title('MITgcm vertical velocity  (m/d)','FontSize',12,'FontWeight','normal');
ylabel('Depth (m)')
xlabel('Latitude (^oN)');
set(gca,'FontSize',12);
colormap(gca,cmocean('balance',20))
hold on;
caxis(caxis);
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(linspace(-8,8,4),-50.*ones(4,1),'k-','linewidth',2)
