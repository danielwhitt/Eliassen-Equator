%% Dan Whitt (daniel.b.whitt@nasa.gov) 
%% Copyright Dan Whitt
%% Written With Matlab v2023b 
%% Upload May 19 2025
%% dependencies
%./data/mitgcm_grid.mat
%./data/mitgcm20yr_buoyavg.nc
%./data/File_y1999_y2018_ub.nc
%./data/wind_stress_zonal_monthly_maps.nc
%./data/RG_ArgoClim_33pfit_2019_mean.nc
%./data/adcp0n110w_mon.nc
%./data/adcp0n140w_mon.nc
%./data/adcp0n170w_mon.nc
%./data/johnson-eq-pac-adcp.nc

clear all;
close all;
id='MATLAB:polyfit:RepeatedPointsOrRescale'
disp('Wyrtki diagram direct simulation numbers...')
vmit3d=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','v'));
umit3d=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','u'));
wmit3d=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','w'));

load('./data/mitgcm_grid.mat','XC','YC','DRF');

XC=XC(41:1460,41:440)-170;
YC=YC(41:1460,41:440);
DXC=zeros(size(XC));
DYC=zeros(size(YC));
DXC(2:end-1,2:end-1)=111300.*cosd(YC(2:end-1,2:end-1)).*(XC(3:end,2:end-1)-XC(1:end-2,2:end-1))./2;
DYC(2:end-1,2:end-1)=111300.*(YC(2:end-1,3:end)-YC(2:end-1,1:end-2))./2;
DXC(1,:)=DXC(2,:);
DXC(end,:)=DXC(end-1,:);
DXC(:,1)=DXC(:,2);
DXC(:,end)=DXC(:,end-1);
DXC(1,1)=DXC(2,2);
DXC(1,end)=DXC(2,end-1);
DXC(end,end)=DXC(end-1,end-1);
DXC(end,1)=DXC(end-1,2);
DYC(1,:)=DYC(2,:);
DYC(end,:)=DYC(end-1,:);
DYC(:,1)=DYC(:,2);
DYC(:,end)=DYC(:,end-1);
DYC(1,1)=DYC(2,2);
DYC(1,end)=DYC(2,end-1);
DYC(end,end)=DYC(end-1,end-1);
DYC(end,1)=DYC(end-1,2);
latmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','latitude'));
lonmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','longitude'));


disp('gcm direct...')

Vmit50m5sdirect=mean(DXC(10,100:101),2).*DRF(1,1,1).*double(sum(sum(mean(vmit3d(:,100:101,1:20),2),1),3))./1e6
Vmit50m5ndirect=mean(DXC(10,300:301),2).*DRF(1,1,1).*double(sum(sum(mean(vmit3d(:,300:301,1:20),2),1),3))./1e6

Umit50m97wdirect=DYC(10,101).*DRF(1,1,1).*double(sum(sum(umit3d(end,101:300,1:20),2),3))./1e6;
Umit50m168wdirect=DYC(10,100).*DRF(1,1,1).*double(sum(sum(umit3d(1,101:300,1:20),2),3))./1e6;
zondiv50m97minus168direct=Umit50m97wdirect-Umit50m168wdirect

hordivmit50mdirect=+Umit50m97wdirect-Umit50m168wdirect+Vmit50m5ndirect-Vmit50m5sdirect;

Wmit50mdirect=double(nansum(nansum(DYC(:,101:300).*DXC(:,101:300).*squeeze(nanmean(wmit3d(:,101:300,21),3)),1),2))./1e6


Vmit50to200m5sdirect=mean(DXC(10,100:101),2).*DRF(1,1,1).*double(sum(sum(mean(vmit3d(:,100:101,21:80),2),1),3))./1e6
Vmit50to200m5ndirect=mean(DXC(10,300:301),2).*DRF(1,1,1).*double(sum(sum(mean(vmit3d(:,300:301,21:80),2),1),3))./1e6

Umit50to200m97wdirect=DYC(10,101).*DRF(1,1,1).*double(sum(sum(umit3d(end,101:300,21:80),2),3))./1e6;
Umit50to200m168wdirect=DYC(10,100).*DRF(1,1,1).*double(sum(sum(umit3d(1,101:300,21:80),2),3))./1e6;
zondiv50to200m97minus168direct=Umit50to200m97wdirect-Umit50to200m168wdirect

hordivmit50to200mdirect=+Umit50to200m97wdirect-Umit50to200m168wdirect+Vmit50to200m5ndirect-Vmit50to200m5sdirect;

Wmit200mdirect=double(nansum(nansum(DYC(:,101:300).*DXC(:,101:300).*squeeze(nanmean(wmit3d(:,101:300,81),3)),1),2))./1e6

verdivmit50to200mdirect=Wmit50mdirect-Wmit200mdirect;


clear vmit3d wmit3d umit3d

disp('gcm mass balance...')
disp('gcm Ekman terms...')

rhoref = 1035; % rho_0
g=9.81; 
tauxmit=double(squeeze(ncread('./data/File_y1999_y2018_ub.nc','Um_Ext')));
tauxmit=squeeze(tauxmit(:,:,1)).*rhoref.*DRF(1,1,1);
f=7.29e-5.*2.*sind(latmit);
fg=repmat(f',[1420 1]);

Ekmantransport5ngcm=sum(DXC(10,300).*mean(-tauxmit(:,300:301)./rhoref./fg(:,300:301),2),1)./1e6
Ekmantransport5sgcm=sum(DXC(10,300).*mean(-tauxmit(:,100:101)./rhoref./fg(:,100:101),2),1)./1e6

disp('gcm mass balance...')
disp('gcm hydrography terms...')

zmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','depth'));
Tmit3d=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','theta')));
Smit3d=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','salt')));

Pmit=sw_pres(abs(zmit),0);

ITmit3d = sw_temp(Smit3d,Tmit3d,repmat(reshape(Pmit,[1 1 136]),[1420 400 1]),zeros(size(Tmit3d)));
Dmit3d = calculate_dynamic_height_seawater(ITmit3d, Smit3d, squeeze(Pmit));
clear ITmit3d

Dmit500mref=-Dmit3d(:,:,1:125)+repmat(squeeze(mean(Dmit3d(:,:,125:126),3)),[1 1 125]);
% g/f int D97-D168 dz
Vmit50m5s=DRF(1,1,1).*g./mean(f(100:101))...
    .*mean(sum(squeeze(mean(Dmit500mref((end-100):end,100:101,1:20),1)-mean(Dmit500mref(1:100,100:101,1:20),1))',1),2)./1e6./1e2.*71/66 % 1e2 converts cm to m in dyn-height, 71/66 rescales for the 5 deg windows on each side
Vmit50m5n=DRF(1,1,1).*g./mean(f(300:301))...
    .*mean(sum(squeeze(mean(Dmit500mref((end-100):end,300:301,1:20),1)-mean(Dmit500mref(1:100,300:301,1:20),1))',1),2)./1e6./1e2.*71/66 % 1e2 converts cm to m in dyn-height, 71/66 rescales for the 5 deg windows on each side

Vmit50to200m5s=DRF(1,1,1).*g./mean(f(100:101))...
        .*mean(sum(squeeze(mean(Dmit500mref((end-100):end,100:101,21:80),1)-mean(Dmit500mref(1:100,100:101,21:80),1))',1),2)./1e6./1e2.*71/66 % 1e2 converts cm to m in dyn-height, 71/66 rescales for the 5 deg windows on each side
Vmit50to200m5n=DRF(1,1,1).*g./mean(f(300:301))...
        .*mean(sum(squeeze(mean(Dmit500mref((end-100):end,300:301,21:80),1)-mean(Dmit500mref(1:100,300:301,21:80),1))',1),2)./1e6./1e2.*71/66 % 1e2 converts cm to m in dyn-height, 71/66 rescales for the 5 deg windows on each side



Wmit50massbalance=Ekmantransport5ngcm-Ekmantransport5sgcm+Vmit50m5n-Vmit50m5s+zondiv50m97minus168direct
hordiv50to200mmassbalance=Vmit50to200m5n-Vmit50to200m5s+zondiv50to200m97minus168direct;
Wmit200mmassbalance=Wmit50massbalance+hordiv50to200mmassbalance

Vmit5n50mmassbalance=Vmit50m5n+Ekmantransport5ngcm
Vmit5s50mmassbalance=Vmit50m5s+Ekmantransport5sgcm

clear Smit3d Tmit3d 

%% observations
clear all;
close all;
id='MATLAB:polyfit:RepeatedPointsOrRescale'

disp('Risien and Chelton Ekman transports...')
addpath ~/'OneDrive - NASA'/software/rainevent-master/cmocean/
monms={'january','february','march','april','may','june','july','august','september','october','november','december'};
latitudewstress=ncread('./data/wind_stress_zonal_monthly_maps.nc','latitude',[241],[80]);
longitudewstress=ncread('./data/wind_stress_zonal_monthly_maps.nc','longitude',[769],[284]);
longitudewstress=longitudewstress-360;

taux=zeros(length(longitudewstress),length(latitudewstress),12);
for mo = 1:12
taux(:,:,mo)=ncread('./data/wind_stress_zonal_monthly_maps.nc',monms{mo},[769 241],[284 80]);
end

[latg,long]=meshgrid(latitudewstress,longitudewstress);
fg=7.29e-5.*2.*sind(latg);

taux(taux==-9999)=nan;
taux=squeeze(mean(taux,3));

DXR=zeros(size(latg));
DXR(2:end-1,2:end-1)=111300.*cosd(latg(2:end-1,2:end-1)).*(long(3:end,2:end-1)-long(1:end-2,2:end-1))./2;
DXR(1,:)=DXR(2,:);
DXR(end,:)=DXR(end-1,:);
DXR(:,1)=DXR(:,2);
DXR(:,end)=DXR(:,end-1);
rhoref=1029;

Ekmantransport5nrisien=sum(mean(DXR(:,60:61).*(-taux(:,60:61))./rhoref./fg(:,60:61),2),1)./1e6
Ekmantransport5srisien=sum(mean(DXR(:,60:61).*(-taux(:,20:21))./rhoref./fg(:,20:21),2),1)./1e6

%
disp('Argo transports...')
g = 9.81; %gravity
rhoref = 1029; % rho_0

LATargo=smooth(double(squeeze(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','LATITUDE',[331],[120])))',5)';
LONargo=ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','LONGITUDE',[1033],[426]);
LONargo=360-LONargo;
LONargo=-LONargo;
Pargo=double(repmat(reshape(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','PRESSURE'),[1 1 58]),[426 120 1]));
Targo=double(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','ARGO_TEMPERATURE_MEAN',[1033 331 1],[426 120 58]));
Sargo=double(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','ARGO_SALINITY_MEAN',[1033 331 1],[426 120 58]));

PDargo=sw_pden(Sargo,Targo,Pargo,zeros(size(Pargo)));
THargo=sw_ptmp(Sargo,Targo,Pargo,zeros(size(Pargo)));
Zargo=-sw_dpth(squeeze(Pargo(1,:,:)),repmat(LATargo,[1 58])); 
Dargo = calculate_dynamic_height_seawater(Targo, Sargo, squeeze(Pargo(1,1,:)));
Dargo500mref=-Dargo(:,:,1:34)+repmat(squeeze(Dargo(:,:,34)),[1 1 34]);

clear latg long
[latg,long]=meshgrid(LATargo,LONargo);
fg=7.29e-5.*2.*sind(latg);

DYC=zeros(size(latg));
DYC(:,2:end-1)=111300.*(latg(:,3:end)-latg(:,1:end-2))./2;
DYC(:,1)=DYC(:,2);
DYC(:,end)=DYC(:,end-1);


% g/f int D97-D168 dz
Vargo50m5s=g./squeeze(mean(fg(20,30:31),2)).*...
    trapz(abs(Zargo(60,1:6)),squeeze(mean(mean(Dargo500mref((end-19):end,30:31,1:6),1),2)-mean(mean(Dargo500mref(1:20,30:31,1:6),1),2))')./1e6./1e2.*71/66+...
    +2.5.*g./squeeze(mean(fg(20,30:31),2)).*...
     (squeeze(mean(mean(Dargo500mref((end-19):end,30:31,1),1),2)-mean(mean(Dargo500mref(1:20,30:31,1),1),2))')./1e6./1e2.*71/66

Vargo50m5n=g./squeeze(mean(fg(20,90:91),2)).*...
    trapz(abs(Zargo(60,1:6)),squeeze(mean(mean(Dargo500mref((end-19):end,90:91,1:6),1),2)-mean(mean(Dargo500mref(1:20,90:91,1:6),1),2))')./1e6./1e2.*71/66+...
    +2.5.*g./squeeze(mean(fg(20,90:91),2)).*...
     (squeeze(mean(mean(Dargo500mref((end-19):end,90:91,1),1),2)-mean(mean(Dargo500mref(1:20,90:91,1),1),2))')./1e6./1e2.*71/66


Vargo50to200m5s=g./squeeze(mean(fg(20,30:31),2)).*...
    trapz(abs(Zargo(60,6:20)),squeeze(mean(mean(Dargo500mref((end-19):end,30:31,6:20),1),2)-mean(mean(Dargo500mref(1:20,30:31,6:20),1),2))')./1e6./1e2.*71/66+...
    +0.9.*g./squeeze(mean(fg(20,30:31),2)).*...
     (squeeze(mean(mean(Dargo500mref((end-19):end,30:31,20),1),2)-mean(mean(Dargo500mref(1:20,30:31,20),1),2))')./1e6./1e2.*71/66

Vargo50to200m5n=g./squeeze(mean(fg(20,90:91),2)).*...
    trapz(abs(Zargo(60,6:20)),squeeze(mean(mean(Dargo500mref((end-19):end,90:91,6:20),1),2)-mean(mean(Dargo500mref(1:20,90:91,6:20),1),2))')./1e6./1e2.*71/66+...
    +0.9.*g./squeeze(mean(fg(20,90:91),2)).*...
     (squeeze(mean(mean(Dargo500mref((end-19):end,90:91,20),1),2)-mean(mean(Dargo500mref(1:20,90:91,20),1),2))')./1e6./1e2.*71/66


% 
disp('zonal convergences from Argo...')
warning('off',id)

for i = 1:120
    for j = 1:34
        clear p
        p=polyfit(squeeze(LONargo-mean(LONargo)),squeeze(Dargo500mref(:,i,j)),3);
        Dargo500mrefp(:,i,j)=polyval(p,LONargo-mean(LONargo));        
    end
end
warning('on',id)



Ugargo = zeros(size(Dargo500mref));
Ugargo(:,2:end-1,:)=-g.*(Dargo500mrefp(:,3:end,:)-Dargo500mrefp(:,1:end-2,:))./(2.*DYC(:,2:end-1))./fg(:,2:end-1)./1e2;
xidx=[43:53,58:63,68:78];
xidx2=[43:53,68:78];

utao110=squeeze(double(ncread('./data/adcp0n110w_mon.nc','u_1205')));
ztao110=ncread('./data/adcp0n110w_mon.nc','depth');
utao110=utao110./100;
utao110(utao110>1e1)=nan;
utao110=nanmean(utao110,2);
utao140=squeeze(double(ncread('./data/adcp0n140w_mon.nc','u_1205')));
ztao140=ncread('./data/adcp0n140w_mon.nc','depth');
utao140=utao140./100;
utao140(utao140>1e1)=nan;
utao140=nanmean(utao140,2);
utao170=squeeze(double(ncread('./data/adcp0n170w_mon.nc','u_1205')));
ztao170=ncread('./data/adcp0n170w_mon.nc','depth');
utao170=utao170./100;
utao170(utao170>1e1)=nan;
utao170=nanmean(utao170,2);
utao170(1:4)=nan;
utao140(1:4)=nan;
utao110(1:4)=nan;
utao170(ztao170>275)=nan;
utao140(ztao140>275)=nan;
utao110(ztao110>275)=nan;
utao110argo=interp1(double(ztao110(~isnan(utao110))),double(utao110(~isnan(utao110))),double(-Zargo(10,1:34)));
utao140argo=interp1(double(ztao140(~isnan(utao140))),double(utao140(~isnan(utao140))),double(-Zargo(10,1:34)));
utao170argo=interp1(double(ztao170(~isnan(utao170))),double(utao170(~isnan(utao170))),double(-Zargo(10,1:34)));

warning('off',id)
for j = 1:34
    clear p
    if ~isnan(utao140argo(j))
    p=polyfit([-170 -140 -110]-mean(LONargo),[utao170argo(j) utao140argo(j) utao110argo(j)],2);
    Ugargo(:,58:63,j)=repmat(polyval(p,LONargo-mean(LONargo)),[1 6]);
    end
end


for i = 1:426
    for j = 1:34
        clear p
        if ~isnan(utao140argo(j))
        p=polyfit(LATargo(xidx),squeeze(Ugargo(i,xidx,j)),6);
        else
        p=polyfit(LATargo(xidx2),squeeze(Ugargo(i,xidx2,j)),3);
        end
        Ugargo(i,53:68,j)=polyval(p,LATargo(53:68));
        polytemp=polyval(p,LATargo(43:52))';
        w=(LATargo(43:52)'+3)./1.75;
        Ugargo(i,43:52,j)=squeeze(Ugargo(i,43:52,j)).*(1-w)+w.*polytemp;
        polytemp=polyval(p,LATargo(69:78))';
        w=(3-LATargo(69:78)')./1.75;
        Ugargo(i,69:78,j)=squeeze(Ugargo(i,69:78,j)).*(1-w)+w.*polytemp;
    end
end


for i = 1:120
    for j = 1:34
        clear p
        p=polyfit(LONargo-mean(LONargo),squeeze(Ugargo(:,i,j)),1);
        dUdxargo(i,j)=p(1)./110750;
        
    end
end
warning('on',id)


figure;
subplot(2,1,1),...
contourf(LATargo',Zargo(10,1:34)',squeeze(mean(Ugargo,1))',linspace(-2,2,41));
colormap(jet(20));
caxis([-1 1]);
ylim([-300 0])
xlim([-8 8])
cbh=colorbar;
xlabel('latitude')
ylabel('depth (m)')
set(gcf,'color','w')
title('(a) Argo zonal mean U')
ylabel(cbh,'ms^{-1} ')
set(gca,'fontsize',14)

subplot(2,1,2),...
contourf(LATargo',Zargo(10,1:34)',(dUdxargo.*1e8)',linspace(-20,20,41));
colormap(jet(20));
caxis([-10 10]);
ylim([-300 0])
xlim([-8 8])
cbh=colorbar;
xlabel('latitude')
ylabel('depth (m)')
set(gcf,'color','w')
title('(b) Argo zonal mean dU/dx')
ylabel(cbh,'10^{-8} s^{-1} ')
set(gca,'fontsize',14)
print(gcf,'U_sections_argo.png','-dpng','-r200')

dUdxargo0to50m=.1666.*111300.*sum(trapz(abs(Zargo(10,1:6)),dUdxargo(31:90,1:6),2),1).*71*110750./1e6+...
2.5.*.1666.*111300.*sum(dUdxargo(31:90,1),1).*71*110750./1e6
dUdxargo50to200m=.1666.*111300.*sum(trapz(abs(Zargo(10,6:20)),dUdxargo(31:90,6:20),2),1).*71*110750./1e6+...
.9.*.1666.*111300.*sum(dUdxargo(31:90,20),1).*71*110750./1e6

%

disp('Zonal convergences from Johnson ADCP data...')
Ujohnson=ncread('./data/johnson-eq-pac-adcp.nc','UM');

PDjohnson=ncread('./data/johnson-eq-pac-adcp.nc','SIGMAM');
XLONjohnson=ncread('./data/johnson-eq-pac-adcp.nc','XLON');
XLONjohnson=XLONjohnson(5:10);
YLATjohnson=ncread('./data/johnson-eq-pac-adcp.nc','YLAT11_101');
dUdxjohnson=zeros(91,50);
warning('off',id)

for i = 1:91
    for j=1:50
        clear p
        p=polyfit(squeeze(XLONjohnson-mean(XLONjohnson)),squeeze(Ujohnson(5:10,i,j)),1);
        dUdxjohnson(i,j)=p(1)./110750;
        clear p
        p=polyfit(XLONjohnson-mean(XLONjohnson),Ujohnson(5:10,i,j),3);
        Upjohnson(i,j)=mean(polyval(p,linspace(XLONjohnson(1)-mean(XLONjohnson),XLONjohnson(end)-mean(XLONjohnson),201)));
    end
end
warning('on',id)

ZDEPjohnson=ncread('./data/johnson-eq-pac-adcp.nc','ZDEP1_50');
Ujohnson=squeeze(mean(Ujohnson(5:10,:,:),1,'omitnan'));
PDjohnson=squeeze(mean(PDjohnson(5:10,:,:),1,'omitnan'));


dUdxadcp0to50m=trapz(YLATjohnson(16:66).*111300,sum(dUdxjohnson(16:66,1:5),2).*10).*71*110750./1e6
dUdxadcp50to200m=trapz(YLATjohnson(16:66).*111300,sum(dUdxjohnson(16:66,6:20),2).*10).*71*110750./1e6


Vobs50m5nmassbalance=Ekmantransport5nrisien+Vargo50m5n
Vobs50m5smassbalance=Ekmantransport5srisien+Vargo50m5s

Wobs50mmassbalance=Vobs50m5nmassbalance-Vobs50m5smassbalance+dUdxadcp0to50m

Wobs200mmassbalance=Wobs50mmassbalance+Vargo50to200m5n-Vargo50to200m5s+dUdxadcp50to200m

Wobs50mmassbalanceargo=Vobs50m5nmassbalance-Vobs50m5smassbalance+dUdxargo0to50m

Wobs200mmassbalanceargo=Wobs50mmassbalanceargo+Vargo50to200m5n-Vargo50to200m5s+dUdxargo50to200m




