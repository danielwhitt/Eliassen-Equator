%% Dan Whitt (daniel.b.whitt@nasa.gov) 
%% Copyright Dan Whitt
%% Written With Matlab v2023b 
%% Github upload May 19 2025
% dependencies
% cmocean 
% ./data/NOAA_GDP310_drifter_monthlymeans.nc
% ./data/File_y1999-2018_buoy_monclim_zonmean.nc
% ./data/sst.ltm.1991-2020.nc
%./data/File_y1999-2018_monclim_taux.nc
% ./data/File_y1999-2018_buoy_monmeans_zonmean.nc
% mitgcm_offline_KE_monthly_zonmean.mat
% ./data/wmitmapmonclimmay19.mat
clear all;
close all;
addpath ./cmocean/
stidxlat=253
lenidxlat=80
stidxlon=49
lenidxlon=284
pathdr='./data/NOAA_GDP310_drifter_monthlymeans.nc'
V=squeeze(double(ncread(pathdr,'V',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])));
U=squeeze(double(ncread(pathdr,'U',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])));

Vdr=squeeze(mean(double(ncread(pathdr,'V',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])),2));

LatDr=squeeze(double(ncread(pathdr,'Lat',[stidxlat],[lenidxlat ])));
LonDr=squeeze(double(ncread(pathdr,'Lon',[stidxlon],[lenidxlon ])));

for it=1:12
    Vdt=detrend(squeeze(V(:,:,it))',2)';
    Vp=V-Vdt;
    stdVdt(:,it)=std(Vdt,1,2);
    for i = 1:80
        [r,l]=xcorr(Vdt(i,:)','normalized');
        rr(:,i)=r;
        [m,lidx]=min(abs(r(284:304)));
        lidx=lidx-1;
        neff(i)=284./lidx;
    end
    neff=neff';
    neff=smooth(neff,5)';
    neff=smooth(neff,5)';
    eVme(:,it)=stdVdt(:,it)./sqrt(neff);
end

for it=1:12
    Udt=detrend(squeeze(U(:,:,it))',2)';
    Up=U-Udt;
    stdUdt(:,it)=std(Udt,1,2);
    for i = 1:80
        [r,l]=xcorr(Udt(i,:)','normalized');
        rr(:,i)=r;
        [m,lidx]=min(abs(r(284:304)));
        lidx=lidx-1;
        neff(i)=284./lidx;
    end
    neff=neff';
    neff=smooth(neff,5)';
    neff=smooth(neff,5)';
    eUme(:,it)=stdUdt(:,it)./sqrt(neff);
end



v15mmit=squeeze(mean(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','v',[1 1 6 1],[1 400 2 12]),3));
v1p25mmit=squeeze(mean(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','v',[1 1 1 1],[1 400 1 12]),3));

u15mmit=squeeze(mean(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','u',[1 1 6 1],[1 400 2 12]),3));

latmit=ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','y',[1],[400]);
zmit=ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','depth',[1],[136]);

w50mmit=86400.*squeeze(mean(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','w',[1 1 20 1],[1 400 2 12]),3));


figure('position',[50 50 800 800]);
subplot(3,2,1),...
contourf(1:12,LatDr,squeeze(mean(V,2))-repmat(squeeze(mean(mean(V,2),3)),[1 12]),-.1:.005:.1);
colormap(gca,cmocean('balance',24));
caxis([-0.06 .06])
cbh=colorbar;
ylabel(cbh,'m/s')
ylim([-8 8])
grid on
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title(['(a) Climatological seasonal cycle      ';' meridional velocity at 15 m (drifters)']);
set(gca,'fontsize',13)

subplot(3,2,3),...
contourf(1:12,latmit,squeeze(v15mmit)-repmat(squeeze(mean(v15mmit,2)),[1 12]),-.1:.005:.1);
colormap(gca,cmocean('balance',24));
caxis([-0.06 .06])
cbh=colorbar;
ylabel(cbh,'m/s')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(b) MITgcm')
set(gca,'fontsize',13)
ylim([-8 8])
grid on
subplot(3,2,5),...
contourf(1:12,LatDr,2.*eVme,linspace(0,.05,10));
caxis([0 .03])
cbh=colorbar;
ylabel(cbh,'m/s')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(c) Obs uncertainty (2 std errors)')
set(gca,'fontsize',13)
set(gcf,'color','w');
ylim([-8 8])
colormap(gca,jet(6));
grid on

subplot(3,2,2),...
contourf(1:12,LatDr,squeeze(mean(U,2))-repmat(squeeze(mean(mean(U,2),3)),[1 12]),-1:.05:1);
colormap(gca,cmocean('balance',20));
caxis([-0.5 .5])
cbh=colorbar;
ylabel(cbh,'m/s')
ylim([-8 8])

xlabel('Month')
ylabel('Latitude ^{\circ}N')
title(['(d) Climatological seasonal cycle  ';'  zonal velocity at 15 m (drifters)']);
set(gca,'fontsize',13)
grid on

subplot(3,2,4),...
contourf(1:12,latmit,squeeze(u15mmit)-repmat(squeeze(mean(u15mmit,2)),[1 12]),-1:.05:1);
colormap(gca,cmocean('balance',20));
caxis([-0.5 .5])
cbh=colorbar;
ylabel(cbh,'m/s')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(e) MITgcm')
set(gca,'fontsize',13)
ylim([-8 8])
grid on


subplot(3,2,6),...
contourf(1:12,LatDr,2.*eUme,linspace(0,.2,20));
colormap(gca,jet(10));
caxis([0 .1])
cbh=colorbar;
ylabel(cbh,'m/s')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(f) Obs uncertainty (2 std errors)')
set(gca,'fontsize',13)
set(gcf,'color','w');
ylim([-8 8])
grid on
print(gcf,'MITgcmclim_seasonalcycle_UandV.png','-dpng','-r200')


for it = 1:12
    pp=polyfit(latmit,v15mmit(:,it),9);
    v15mmitp(:,it)=polyval(pp,latmit);
    pp=polyfit(LatDr,squeeze(mean(V(:,:,it),2)),9);
    Vdrp(:,it)=polyval(pp,LatDr);
end




dVdymit=zeros(size(v15mmit));
dVdymit1p25=zeros(size(v15mmit));
dVdymit(2:end-1,:)=86400.*(v15mmit(3:end,:)-v15mmit(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));
dVdymit1p25(2:end-1,:)=86400.*(v1p25mmit(3:end,:)-v1p25mmit(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));



figure;

subplot(1,2,1),...
contourf(1:12,latmit,squeeze(dVdymit),-.04:.005:.04);
caxis([-.04 .04])
cbh=colorbar;
ylabel(cbh,'1/d')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(b) dV/dy at 15 m')
set(gca,'fontsize',13)
ylim([-6 6])
grid on
set(gcf,'color','w')
[m,latmidx]=max(dVdymit,[],1);
hold on
plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,cmocean('balance',16));

subplot(1,2,2),...
contourf(1:12,latmit,squeeze(w50mmit),-5:.2:5);
caxis([-1.6 1.6])
cbh=colorbar;
ylabel(cbh,'m/d')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(a) w at 50 m')
set(gca,'fontsize',13)
ylim([-6 6])
grid on
[m,latmidx]=max(w50mmit,[],1);
hold on
plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,cmocean('balance',16));
print(gcf,'MITgcmclim_seasonalcycle_dVdyW.png','-dpng','-r200')



figure;

subplot(1,2,1),...
contourf(1:12,latmit,squeeze(dVdymit1p25),-.04:.005:.04);
caxis([-.04 .04])
cbh=colorbar;
ylabel(cbh,'1/d')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(b) dV/dy at 1.25 m')
set(gca,'fontsize',13)
ylim([-6 6])
grid on
set(gcf,'color','w')
[m,latmidx]=max(dVdymit1p25,[],1);
hold on
plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,cmocean('balance',16));

subplot(1,2,2),...
contourf(1:12,latmit,squeeze(w50mmit),-5:.2:5);
caxis([-1.6 1.6])
cbh=colorbar;
ylabel(cbh,'m/d')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(a) w at 50 m')
set(gca,'fontsize',13)
ylim([-6 6])
grid on
[m,latmidx]=max(w50mmit,[],1);
hold on
plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,cmocean('balance',16));
print(gcf,'MITgcmclim_seasonalcycle_dVdy1p25W.png','-dpng','-r200')



w50mmitsa=w50mmit-repmat(mean(w50mmit,2),[1 12]);
dVdymitsa=dVdymit-repmat(mean(dVdymit,2),[1 12]);
w50mmitsa=w50mmitsa(121:320,:);
dVdymitsa=dVdymitsa(121:320,:);


figure;
subplot(2,1,2),...
scatter(dVdymitsa(:),w50mmitsa(:),10,'k.');
[r,p]=corrcoef(dVdymitsa(:),w50mmitsa(:));
rr=r(1,2);
xlabel('dV/dy at 15 m (1/d)')
ylabel('w at 50 m (m/d)')
[b,bi,r,ri,stats]=regress(w50mmitsa(:) ,[ones(size(dVdymitsa(:))) dVdymitsa(:) ]);
title(['(b) spatiotemporal correlation, seasonal dV/dy(15 m) vs. w(50 m), r^2=',num2str(stats(1),2),' s=',num2str(b(2),2),' m'],'fontweight','normal');
set(gca,'fontsize',12)
grid on;
set(gcf,'color','w')

dVdymitc=mean(dVdymit,2);
dVdymitc=dVdymitc(121:320);
w50mmitc=mean(w50mmit,2);
w50mmitc=w50mmitc(121:320);


subplot(2,1,1),...
scatter(dVdymitc(:),w50mmitc(:),10,'k.');
xlabel('dV/dy at 15 m (1/d)')
ylabel('w at 50 m (m/d)')
[b,bi,r,ri,stats]=regress(w50mmitc(:) ,[ones(size(dVdymitc(:))) dVdymitc(:) ]);
title(['(a) spatial correlation (4S-6N), mean dV/dy(15 m) vs. w(50 m), r^2=',num2str(stats(1),2),' s=',num2str(b(2),2),' m'],'fontweight','normal');
set(gca,'fontsize',12)
grid on;
set(gcf,'color','w')

print(gcf,'MITgcmclim_seasonalcycle_dVdyW_scatter.png','-dpng','-r200')






SSTmit=squeeze(mean(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','theta',[1 1 1 1],[1 400 1 12]),3));
figure;
subplot(1,2,1),...
contourf(1:12,latmit,squeeze(SSTmit),20:.25:30);
colormap(jet(40));
caxis([20 30])
cbh=colorbar;
ylabel(cbh,'SST deg C')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('MITgcm')
set(gca,'fontsize',13)
ylim([-8 8])
grid on


latsst=ncread('./data/sst.ltm.1991-2020.nc','lat',[81],[20]);
lonsst=ncread('./data/sst.ltm.1991-2020.nc','lon',[193],[71]);
sstltm=squeeze(mean(double(ncread('./data/sst.ltm.1991-2020.nc','sst',[193 81 1],[71 20 12])),1));

subplot(1,2,2),...
contourf(1:12,latsst,squeeze(sstltm),20:.25:30);
colormap(jet(40));
caxis([20 30])
cbh=colorbar;
ylabel(cbh,'SST deg C')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('Reynolds LTM 1990-2020')
set(gca,'fontsize',13)
ylim([-8 8])
grid on


load('./data/mitgcm_offline_KE_monthly_zonmean.mat')

VVmit15m=squeeze(mean(mean(SFnow(:,6:7,:,2,:),5),2));
UUmit15m=squeeze(mean(mean(SFnow(:,6:7,:,1,:),5),2));
UVmit15m=squeeze(mean(mean(SFnow(:,6:7,:,4,:),5),2));

umit=ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','u');
umit15mall=squeeze(mean(umit(1,:,6:7,:),3));
umit15mc=reshape(repmat(mean(reshape(umit15mall,[400 12 20]),3),[1 1 20]),[400 240]);
umit15mp=umit15mall-umit15mc;
timemitmat=double(ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','time')./24+datenum(1950,1,1));

vmit=ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','v');
vmit15mall=squeeze(mean(vmit(1,:,6:7,:),3));
vmit15mc=reshape(repmat(mean(reshape(vmit15mall,[400 12 20]),3),[1 1 20]),[400 240]);
vmit15mp=vmit15mall-vmit15mc;
dVdymit15mp=zeros(size(vmit15mp));
dVdymit15mc=zeros(size(vmit15mc));

Tmit=ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','theta');
SSTmitall=squeeze(mean(Tmit(1,:,1,:),3));
SSTmitc=reshape(repmat(mean(reshape(SSTmitall,[400 12 20]),3),[1 1 20]),[400 240]);
SSTmitp=SSTmitall-SSTmitc;


figure('position',[10 10 1200 300]);
subplot(1,3,1)
contourf(1:12,latmit,squeeze(UUmit15m)-umit15mc(:,1:12).^2,0:.005:.5);
caxis([0 .2])
cbh=colorbar;
ylabel(cbh,'m^2/s^2')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(a) <U^2>-<U>^2 at 15 m')
set(gca,'fontsize',13)
ylim([-8 8])
grid on
set(gcf,'color','w')
%[m,latmidx]=max(dVdymit,[],1);
%hold on
%plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,jet(40));
hold on;

caxis(caxis)
[c,h]=contour(1:12,latmit,squeeze(u15mmit),-1:.1:1,'k-','linewidth',1);
clabel(c,h,'fontsize',12)
hold on
plot(1:12,latmit(latmidx),'r.','markersize',12,'linestyle','none')

subplot(1,3,2)
contourf(1:12,latmit,squeeze(VVmit15m)-vmit15mc(:,1:12).^2,0:.005:.1);
caxis([0 .1])
cbh=colorbar;
ylabel(cbh,'m^2/s^2')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(b) <V^2>-<V>^2 at 15 m')
set(gca,'fontsize',13)
ylim([-8 8])
grid on
set(gcf,'color','w')
%[m,latmidx]=max(dVdymit,[],1);
%hold on
%plot(1:12,latmit(latmidx),'linewidth',2,'color','r')
colormap(gca,jet(20));
hold on;
caxis(caxis)
[c,h]=contour(1:12,latmit,squeeze(u15mmit),-1:.1:1,'k-','linewidth',1);
clabel(c,h,'fontsize',12)
hold on
plot(1:12,latmit(latmidx),'r.','markersize',12,'linestyle','none')


subplot(1,3,3)
contourf(1:12,latmit,squeeze(UVmit15m)-umit15mc(:,1:12).*vmit15mc(:,1:12),-.05:.005:.05);
caxis([-.05 .05])
cbh=colorbar;
ylabel(cbh,'m^2/s^2')
xlabel('Month')
ylabel('Latitude ^{\circ}N')
title('(c) <UV>-<U><V> at 15 m')
set(gca,'fontsize',13)
ylim([-8 8])
grid on
set(gcf,'color','w');
hold on;
caxis(caxis)
[c,h]=contour(1:12,latmit,squeeze(u15mmit),-1:.1:1,'k-','linewidth',1);
clabel(c,h,'fontsize',12)
hold on
plot(1:12,latmit(latmidx),'r.','markersize',12,'linestyle','none')
colormap(gca,cmocean('balance',20));

print(gcf,'MITgcmclim_seasonalcycle_U2andV2andUV.png','-dpng','-r200');



wmit=86400.*ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','w');
Tmit=squeeze(double(ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','theta')));

wmit50mall=squeeze(mean(wmit(1,:,20:21,:),3));
wmit50mc=reshape(repmat(mean(reshape(wmit50mall,[400 12 20]),3),[1 1 20]),[400 240]);
wmit50mp=wmit50mall-wmit50mc;

vmit=ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','v');
vmit15mall=squeeze(mean(vmit(1,:,6:7,:),3));
vmit15mc=reshape(repmat(mean(reshape(vmit15mall,[400 12 20]),3),[1 1 20]),[400 240]);
vmit15mp=vmit15mall-vmit15mc;
dVdymit15mp=zeros(size(vmit15mp));
dVdymit15mc=zeros(size(vmit15mc));

dVdymit15mp(2:end-1,:)=86400.*(vmit15mp(3:end,:)-vmit15mp(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));
dVdymit15mc(2:end-1,:)=86400.*(vmit15mc(3:end,:)-vmit15mc(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));


vmit1p25mall=squeeze(mean(vmit(1,:,1,:),3));
vmit1p25mc=reshape(repmat(mean(reshape(vmit1p25mall,[400 12 20]),3),[1 1 20]),[400 240]);
vmit1p25mp=vmit1p25mall-vmit1p25mc;
dVdymit1p25mp=zeros(size(vmit1p25mp));
dVdymit1p25mp(2:end-1,:)=86400.*(vmit1p25mp(3:end,:)-vmit1p25mp(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));
dVdymit1p25mc=zeros(size(vmit1p25mc));
dVdymit1p25mc(2:end-1,:)=86400.*(vmit1p25mc(3:end,:)-vmit1p25mc(1:end-2,:))./(111000.*(latmit(3)-latmit(1)));



umit=ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','u');
umit15mall=squeeze(mean(umit(1,:,6:7,:),3));
umit15mc=reshape(repmat(mean(reshape(umit15mall,[400 12 20]),3),[1 1 20]),[400 240]);
umit15mp=umit15mall-umit15mc;
timemitmat=double(ncread('./data/File_y1999-2018_buoy_monmeans_zonmean.nc','time')./24+datenum(1950,1,1));




figure('Position',[50 50 1200 400]);
subplot(2,3,[1 2]),...
contourf(timemitmat,latmit,wmit50mp,-3:.25:3,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-1.5 1.5])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(a) Aseasonal monthly vertical velocity anomaly at 50 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',24))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),1.5.*ones(3,1),'magenta-','linewidth',1)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),2.5.*ones(3,1),'magenta-','linewidth',1)
cbh=colorbar;
ylabel(cbh,'m/d')
grid on



subplot(2,3,[4 5]),...
contourf(timemitmat,latmit,dVdymit15mp,-.05:.005:.05,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-.05 .05])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(b) Aseasonal dV/dy at 15 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',20))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
cbh=colorbar;
ylabel(cbh,'1/d')
grid on

wmit50mpres=wmit50mp(121:320,:);
dVdymit15mpres=dVdymit15mp(121:320,:);
[r,p]=corrcoef(wmit50mpres(:),dVdymit15mpres(:));
rrdVdy=r(1,2);
subplot(2,3,[3 6]),...
scatter(dVdymit15mpres(:),wmit50mpres(:),10,'k','filled');
%datetick(cbh)
xlabel('Aseasonal monthly dV/dy at 15 m (1/d)')
ylabel('Aseasonal monthly w at 50 m (m/d)')
title(['(c) Correlation in time and space from 4^oS-6^oN,r^2=',num2str(rrdVdy.^2,2)])
grid on
axis square
xlim([-.05 .05])
ylim([-1.5 1.5])
set(gcf,'color','w')
set(gca,'fontsize',12)
grid on
print(gcf,'MITgcm_aseasonal_WanddVdy.png','-dpng','-r400')



figure('Position',[50 50 1200 400]);
subplot(2,3,[1 2]),...
contourf(timemitmat,latmit,wmit50mp,-3:.25:3,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-1.5 1.5])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(a) Aseasonal monthly vertical velocity anomaly at 50 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',24))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),1.5.*ones(3,1),'magenta-','linewidth',1)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),2.5.*ones(3,1),'magenta-','linewidth',1)
cbh=colorbar;
ylabel(cbh,'m/d')
grid on



subplot(2,3,[4 5]),...
contourf(timemitmat,latmit,dVdymit1p25mp,-.05:.005:.05,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-.05 .05])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(b) Aseasonal dV/dy at 1.25 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',20))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
cbh=colorbar;
ylabel(cbh,'1/d')
grid on

wmit50mpres=wmit50mp(121:320,:);
dVdymit1p25mpres=dVdymit1p25mp(121:320,:);
[r,p]=corrcoef(wmit50mpres(:),dVdymit1p25mpres(:));
rrdVdy=r(1,2);
subplot(2,3,[3 6]),...
scatter(dVdymit1p25mpres(:),wmit50mpres(:),10,'k','filled');
%datetick(cbh)
xlabel('Aseasonal monthly dV/dy at 1.25 m (top grid cell) (1/d)')
ylabel('Aseasonal monthly w at 50 m (m/d)')
title(['(c) Correlation in time and space from 4^oS-6^oN,r^2=',num2str(rrdVdy.^2,2)])
grid on
axis square
xlim([-.05 .05])
ylim([-1.5 1.5])
set(gcf,'color','w')
set(gca,'fontsize',12)
grid on
print(gcf,'MITgcm_aseasonal_WanddVdy1p25m.png','-dpng','-r400')


figure('Position',[50 50 1200 400]);
subplot(2,3,[1 2]),...
contourf(timemitmat,latmit,dVdymit15mp,-.05:.005:.05,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-.05 .05])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(a) Aseasonal monthly dV/dy anomaly at 15 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',24))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),1.5.*ones(3,1),'magenta-','linewidth',1)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),2.5.*ones(3,1),'magenta-','linewidth',1)
cbh=colorbar;
ylabel(cbh,'1/d')
grid on



subplot(2,3,[4 5]),...
contourf(timemitmat,latmit,dVdymit1p25mp,-.05:.005:.05,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-.05 .05])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(b) Aseasonal dV/dy at 1.25 m'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',20))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
cbh=colorbar;
ylabel(cbh,'1/d')
grid on

dVdymit15mpres=dVdymit15mp(121:320,:);
dVdymit1p25mpres=dVdymit1p25mp(121:320,:);
[r,p]=corrcoef(dVdymit15mpres(:),dVdymit1p25mpres(:));
rrdVdy=r(1,2);
subplot(2,3,[3 6]),...
scatter(dVdymit1p25mpres(:),dVdymit15mpres(:),10,'k','filled');
%datetick(cbh)
xlabel('Aseasonal monthly dV/dy at 1.25 m (top grid cell) (1/d)')
ylabel('Aseasonal monthly dV/dy at 15 m (1/d)')
title(['(c) Correlation in time and space from 4^oS-6^oN,r^2=',num2str(rrdVdy.^2,2)])
grid on
axis square
xlim([-.05 .05])
ylim([-.05 .05])
set(gcf,'color','w')
set(gca,'fontsize',12)
grid on
print(gcf,'MITgcm_aseasonal_dVdy15manddVdy1p25m.png','-dpng','-r400')


for j=1:400
    wmit50mpsm(j,:)=smooth(wmit50mp(j,:),9);
    wmit50mallsm(j,:)=wmit50mc(j,:)+wmit50mpsm(j,:);
end


figure('Position',[50 50 1200 900]);
subplot(5,1,[2 3]),...
contourf(timemitmat,latmit,wmit50mall,-5:.25:5,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-2 2])
xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(b) Monthly and zonal mean w_{50}'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',40))
caxis(caxis);
hold on;
hold on
plot(timemitmat,2.*ones(size(timemitmat)),'k-','linewidth',1);
hold on
plot(timemitmat,-2.*ones(size(timemitmat)),'k-','linewidth',1)

[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
clabel(c,h)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),1.5.*ones(3,1),'magenta-','linewidth',1)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),2.5.*ones(3,1),'magenta-','linewidth',1)
cbh=colorbar;
ylabel(cbh,'m/d')
grid on
set(gcf,'color','w')

w50mmitallasym2=zeros(400,240);
w50mmitallsym=zeros(400,240);
w50mmitallasym=zeros(400,240);

w50mmitallsym(1:200,:)=wmit50mall(1:200,:).*0.5+wmit50mall(400:-1:201,:).*.5;
w50mmitallsym(201:400,:)=wmit50mall(200:-1:1,:).*0.5+wmit50mall(201:400,:).*.5;
w50mmitallasym=wmit50mall-w50mmitallsym;
w50mmitallasym2(201:400,:)=w50mmitallsym(201:400,:)-w50mmitallasym(200:-1:1,:);


[~,latmidxa]=max(wmit50mall,[],1);
[~,idxmax2a]=max(w50mmitallasym2(201:400,:),[],1);
lattemp=latmit(201:400);
hold on
plot(timemitmat,latmit(latmidxa),'r.','markersize',12);
%plot(timemitmat,lattemp(idxmax2a),'bo','markersize',12);
%plot(timemitmat,-lattemp(idxmax2a),'b.','markersize',12);


subplot(5,1,[4 5]),...
contourf(timemitmat,latmit,wmit50mp,-3:.25:3,'linestyle','none');
datetick('x')
shading flat
ylim([-6 6])
caxis([-2 2])
hold on
plot(timemitmat,2.*ones(size(timemitmat)),'k-','linewidth',1);
hold on
plot(timemitmat,-2.*ones(size(timemitmat)),'k-','linewidth',1)

xlim([datenum(1999,1,1) datenum(2019,1,1)])
title(['(c) Interannual anomalies in zonal mean w_{50}'])
ylabel('Latitude ^oN')
set(gca,'fontsize',12)
colormap(gca,cmocean('balance',40))
caxis(caxis);
hold on;
[c,h]=contour(timemitmat,latmit,wmit50mc,-1.5:.5:1.5,'k-');
plot(timemitmat,latmit(latmidxa),'r.','markersize',12);

clabel(c,h)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),1.5.*ones(3,1),'magenta-','linewidth',1)
%plot(linspace(datenum(1999,1,1),datenum(2019,1,1),3),2.5.*ones(3,1),'magenta-','linewidth',1)
cbh=colorbar;
ylabel(cbh,'m/d')
grid on
set(gcf,'color','w')





W50p4S4Nmsmooth=smooth(mean(wmit50mp(101:300,:),1)',9);

%W50p4S4Nmsmooth=smooth((sum(repmat(cosd(latmit(101:300)),[1 240]).*wmit50mp(101:300,:),1)./sum(repmat(cosd(latmit(101:300)),[1 240]),1))',9);

SSTp4S4Nsmooth=smooth(mean(SSTmitp(101:300,:),1)',9);

subplot(5,1,1,'Position',[0.13 0.800677966101695 0.724166666666667 0.124322033898305]),...
    yyaxis left
plot(timemitmat,W50p4S4Nmsmooth.*111300.*10.*71.*110500./1e6./86400,'linewidth',2);
hold on

hold on

    datetick('x','keeplimits')
    ylim([-12 12])
    ylabel('Upwelling at 50m (Sv)')
    yyaxis right
plot(timemitmat,SSTp4S4Nsmooth,'linewidth',2);
    datetick('x')
    title('(a) Regionally integrated monthly upwelling at 50 m and SST anomalies (9-mo smoothed)')
    grid on
    ylabel('SST anomaly (deg C)')
    set(gca,'Fontsize',12)

print(gcf,'MITgcm_aseasonal_WandWa.png','-dpng','-r400')

figure;
plot(latmit,mean(wmit50mall(:,198:206),2),latmit,mean(wmit50mall(:,210:218),2),'--','linewidth',2)
hold on
wasym2015=mean(wmit50mall(:,198:206),2);
wasym2016=mean(wmit50mall(:,210:218),2);

wasym2015(201:400)=wasym2015(201:400)-wasym2015(200:-1:1);
wasym2016(201:400)=wasym2016(201:400)-wasym2016(200:-1:1);

plot(latmit(201:400),wasym2015(201:400),'k-',latmit(201:400),wasym2016(201:400),'k--','linewidth',1)
%area(latmit(201:400),[wasym2016(201:400) wasym2015(201:400)])
%plot(latmit(200:-1:1),-wasym2015(201:400),'k-',latmit(200:-1:1),-wasym2016(201:400),'k--','linewidth',1)


legend('June 2015-Feb 2016','June 2016-Feb 2017','Asymmetry 2015','Asymmetry 2016','location','northwest')
grid on
xlim([-6 6])
xlabel('Latitude')
ylabel('w_{50} (m/d)')
title('Zonal mean vertical velocity w_{50} ')
set(gca,'Fontsize',15)



figure;
plot(latmit,35.*mean(dVdymit15mp(:,198:206)+dVdymit15mc(:,198:206),2),latmit,35.*mean(dVdymit15mp(:,210:218)+dVdymit15mc(:,210:218),2),'--','linewidth',2)
hold on
wasym2015=35.*mean(dVdymit15mp(:,198:206)+dVdymit15mc(:,198:206),2);
wasym2016=35.*mean(dVdymit15mp(:,210:218)+dVdymit15mc(:,210:218),2);

wasym2015(201:400)=wasym2015(201:400)-wasym2015(200:-1:1);
wasym2016(201:400)=wasym2016(201:400)-wasym2016(200:-1:1);

plot(latmit(201:400),wasym2015(201:400),'k-',latmit(201:400),wasym2016(201:400),'k--','linewidth',1)
%area(latmit(201:400),[wasym2016(201:400) wasym2015(201:400)])
%plot(latmit(200:-1:1),-wasym2015(201:400),'k-',latmit(200:-1:1),-wasym2016(201:400),'k--','linewidth',1)


legend('June 2015-Feb 2016','June 2016-Feb 2017','Asymmetry 2015','Asymmetry 2016','location','northwest')
grid on
xlim([-6 6])
xlabel('Latitude')
ylabel('w_{50} (m/d)')
title('Zonal mean scaled dVdy at 15 m')
set(gca,'Fontsize',15)


figure;
plot(latmit,35.*mean(dVdymit1p25mp(:,198:206)+dVdymit1p25mc(:,198:206),2),latmit,35.*mean(dVdymit1p25mp(:,210:218)+dVdymit1p25mc(:,210:218),2),'--','linewidth',2)
hold on
wasym2015=35.*mean(dVdymit1p25mp(:,198:206)+dVdymit1p25mc(:,198:206),2);
wasym2016=35.*mean(dVdymit1p25mp(:,210:218)+dVdymit1p25mc(:,210:218),2);

wasym2015(201:400)=wasym2015(201:400)-wasym2015(200:-1:1);
wasym2016(201:400)=wasym2016(201:400)-wasym2016(200:-1:1);

plot(latmit(201:400),wasym2015(201:400),'k-',latmit(201:400),wasym2016(201:400),'k--','linewidth',1)
%area(latmit(201:400),[wasym2016(201:400) wasym2015(201:400)])
%plot(latmit(200:-1:1),-wasym2015(201:400),'k-',latmit(200:-1:1),-wasym2016(201:400),'k--','linewidth',1)


legend('June 2015-Feb 2016','June 2016-Feb 2017','Asymmetry 2015','Asymmetry 2016','location','northwest')
grid on
xlim([-6 6])
xlabel('Latitude')
ylabel('w_{50} (m/d)')
title('Zonal mean scaled dVdy at 1.25 m')
set(gca,'Fontsize',15)





%%


tauxmit=squeeze(mean(double(ncread('./data/File_y1999-2018_monclim_taux.nc','Um_Ext').*2.5.*1029),1));
w50mmitasym2=zeros(400,12);
w50mmitsym=zeros(400,12);
w50mmitasym=zeros(400,12);

w50mmitsym(1:200,:)=w50mmit(1:200,:).*0.5+w50mmit(400:-1:201,:).*.5;
w50mmitsym(201:400,:)=w50mmit(200:-1:1,:).*0.5+w50mmit(201:400,:).*.5;
w50mmitasym=w50mmit-w50mmitsym;
w50mmitasym2(201:400,:)=w50mmitasym(201:400,:)-w50mmitasym(200:-1:1,:);



figure('position',[10 10 1200 350]);

subplot(1,3,1),...
contourf(1:12,latmit,squeeze(w50mmit),-5:.2:5,'linestyle','none');
hold on;
caxis([-1.6 1.6])
cbh=colorbar;
%ylabel(cbh,'m/d')
xlabel('Month')
ylabel('latitude (^{\circ}N)')
title('(a) w at 50 m (m/d)')
set(gca,'fontsize',12)
ylim([-6 6])
grid on
[m,latmidx]=max(w50mmit,[],1);
[~,idxmax2]=max(w50mmitasym2(201:400,:),[],1);
lattemp=latmit(201:400);
hold on
plot(1:12,latmit(latmidx),'r.','markersize',12);
plot(1:12,lattemp(idxmax2),'b.','markersize',12);
plot(1:12,-lattemp(idxmax2),'b.','markersize',12);

caxis(caxis);
[c,h]=contour(1:12,latmit,tauxmit,linspace(-.1,.1,21),'w');
clabel(c,h,'color','w')

caxis(caxis);
[c,h]=contour(1:12,latmit,tauxmit,linspace(-.1,.1,21),'w','linewidth',1);
clabel(c,h,'color','w')
[c,h]=contour(1:12,latmit,umit15mc(:,1:12),linspace(-1,1,21),'k','linewidth',1);
clabel(c,h,'color','k')


colormap(gca,cmocean('balance',16));

subplot(1,3,2),...
contourf(1:12,latmit,squeeze(35.*dVdymit15mc(:,1:12)),-5:.2:5,'linestyle','none');
hold on;
caxis([-1.6 1.6])
cbh=colorbar;
%ylabel(cbh,'m/d')
xlabel('Month')
ylabel('latitude (^{\circ}N)')
title('(b) Scaled meridional divergence dV/dy at 15 m (m/d)')
set(gca,'fontsize',12)
ylim([-6 6])
grid on
%[m,latmidx]=max(w50mmit,[],1);
%[~,idxmax2]=max(w50mmitasym2(201:400,:),[],1);
%lattemp=latmit(201:400);
hold on
plot(1:12,latmit(latmidx),'r.','markersize',12);
plot(1:12,lattemp(idxmax2),'b.','markersize',12);
plot(1:12,-lattemp(idxmax2),'b.','markersize',12);

caxis(caxis);
[c,h]=contour(1:12,latmit,tauxmit,linspace(-.1,.1,21),'w');
clabel(c,h,'color','w')

caxis(caxis);
[c,h]=contour(1:12,latmit,tauxmit,linspace(-.1,.1,21),'w','linewidth',1);
clabel(c,h,'color','w')
[c,h]=contour(1:12,latmit,umit15mc(:,1:12),linspace(-1,1,21),'k','linewidth',1);
clabel(c,h,'color','k')


colormap(gca,cmocean('balance',16));

for j=1:400
    wmit50mpsm(j,:)=smooth(wmit50mp(j,:),9);
end


subplot(1,3,3),...
plot(std(wmit50mp-wmit50mpsm,1,2),latmit,'linewidth',1,'color','k');
hold on;
plot(std(wmit50mpsm,1,2),latmit,'linewidth',1,'color','r');
plot(std(wmit50mc,1,2),latmit,'linewidth',1,'color','b');
%plot(std(wmit50mp-wmit50mpsm,1,2).^2+std(wmit50mpsm,1,2).^2+std(wmit50mc,1,2),latmit,'g',std(wmit50mp+wmit50mc,1,2),latmit,'cyan')

%plot(mean(wmit50mc,2),latmit,'linewidth',1,'color','r');
ylim([-6 6])
grid on
%ylabel('latitude (^{\circ}N')
xlabel('(m/d)')
title('(c) stdevs of monthly anomalies')
set(gca,'fontsize',12)
xlim([0 0.4])
set(gcf,'color','w')
legend('intraseasonal std','interannual std','mean seasonal cycle std','location','southeast')

print(gcf,'MITgcm_seasonalcycleandstdevanoms_WwithU_3pan.png','-dpng','-r200')
%%
wmitmonclim=86400.*squeeze(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','w',[1 1 1 1],[1 400 136 12]));
Tmitmonclim=squeeze(ncread('./data/File_y1999-2018_buoy_monclim_zonmean.nc','theta',[1 1 1 1],[1 400 136 12]));

figure('position',[50 50 1200 900]);
for i = 1:4
subplot(4,2,i),...
contourf(latmit,zmit(2:end).*0.5+zmit(1:end-1).*0.5,squeeze(mean(wmitmonclim(:,2:end,3*(i-1)+(1:3)),3))',linspace(-2,2,21)); shading flat
hold on
caxis([-2 2])
caxis(caxis)
hold on
[c,h]=contour(latmit,zmit,squeeze(mean(Tmitmonclim(:,:,3*(i-1)+(1:3)),3))',-2:2:30,'linewidth',1,'color','k');
clabel(c,h);
ylabel('depth (m)','FontSize',14,'FontWeight','normal');
xlabel('latitude ({\circ}N)','FontSize',14,'FontWeight','normal');
cbh = colorbar();
set(gca,'FontSize',14,'FontWeight','normal');
if i==1
title('(a) JFM zonal mean w (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 2
title('(b) AMJ zonal mean w (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 3
title('(c) JAS zonal mean w (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 4
title('(d) OND zonal mean w (m/d)','fontsize',14,'FontWeight','normal');
end
ylim([-300 0])
hold on;
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
xlim([-6 6])
plot(linspace(-8,8,17),-50.*ones(17,1),'k-','linewidth',2)
grid on
set(gca,'fontsize',14)
colormap(gca,cmocean('balance',20));
end
set(gcf,'color','w')
%print(gcf,'MITgcmclim_seasonalcycle_Wzonmeansection_4panel.png','-dpng','-r200')

%%
load('./data/wmitmapmonclimmay19.mat');

for i = 1:4
    clear wmitmap
    wmitmap=squeeze(mean(wmitmapmonclim(:,:,3*(i-1)+(1:3)),3));
    wmitmap(isnan(wmitmap))=0;
    %for j=1:400
    %    wmitmap(:,j)=smooth(wmitmap(:,j),101);
    %end
subplot(4,2,4+i),contourf(lonmit,latmit,wmitmap',linspace(-4,4,41),'linestyle','none');
caxis([-2 2]);
colormap(gca,cmocean('balance',20));
hold on
plot(lonmit(:),zeros(1420,1),'k-','linewidth',2)
plot(lonmit(:),2.*ones(1420,1),'k-','linewidth',2)
plot(lonmit(:),-2.*ones(1420,1),'k-','linewidth',2)
grid on;
ylim([-6 6])
cbh=colorbar;
ylabel(cbh,'m/d')
if i==1
title('(e) JFM mean w at 50 m (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 2
title('(f) AMJ mean w at 50 m (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 3
title('(g) JAS mean w at 50 m (m/d)','fontsize',14,'FontWeight','normal');
elseif i == 4
title('(h) OND mean w at 50 m (m/d)','fontsize',14,'FontWeight','normal');
end
set(gca,'color','k','linewidth',1)
ylabel('latitude ({\circ}N)','FontSize',14,'FontWeight','normal');
xlabel('longitude ({\circ}E)','FontSize',14,'FontWeight','normal');
set(gca,'fontsize',14)
hold on
clear wmitmapsym wmitmapasym wmitmapasym2
wmitmapsym(:,1:200)=wmitmap(:,1:200).*0.5+wmitmap(:,400:-1:201).*.5;
wmitmapsym(:,201:400)=wmitmap(:,200:-1:1).*0.5+wmitmap(:,201:400).*.5;
wmitmapasym=wmitmap-wmitmapsym;
wmitmapasym2(:,201:400)=wmitmapasym(:,201:400)-wmitmapasym(:,200:-1:1);
clear lattemp idxmax idxmax2 
[~,idxmax]=max(wmitmap,[],2);
[~,idxmax2]=max(wmitmapasym2(:,201:400),[],2);
lattemp=latmit(201:400);
plot(lonmit(1:5:end),latmit(idxmax(1:5:end)),'r.','markersize',5)
plot(lonmit(1:5:end),lattemp(idxmax2(1:5:end)),'b.','markersize',5)
plot(lonmit(1:5:end),-lattemp(idxmax2(1:5:end)),'b.','markersize',5)

% Create arrow
end
set(gcf,'color','w')

%print(gcf,'MITgcmclim_seasonalcycle_W50maps_4panel.png','-dpng','-r200')



