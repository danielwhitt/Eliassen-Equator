%% Dan Whitt (daniel.b.whitt@nasa.gov) 
%% Copyright Dan Whitt
%% Written With Matlab v2023b 
%% Github upload May 19 2025
%% figs 3 and A1-A5
% dependencies
%./data/mitgcm_grid.mat
%./data/mitgcm_offline_KE_avg.mat
%./data/File_y1999_y2018_ub.nc
%./data/File_y1999_y2018_vb.nc
%./data/NOAA_GDP310_drifter_monthlymeans.nc
%./data/johnson-eq-pac-adcp.nc
%./data/adcp0n170w_mon.nc
%./data/adcp0n140w_mon.nc
%./data/adcp0n110w_mon.nc
%./data/RG_ArgoClim_33pfit_2019_mean.nc
%./seawater/
%./cmocean/
clear all;
close all;
addpath ./cmocean/

% model
load('./data/mitgcm_grid.mat','XC','YC','DRF');



latmit=mean(YC,1,'omitnan');
lonmit=mean(XC,2,'omitnan');
lonmit=-170+lonmit(41:1460);
latmit=latmit(41:440);




%%

stidxlat=253
lenidxlat=80
stidxlon=49
lenidxlon=284
pathdr='./data/NOAA_GDP310_drifter_monthlymeans.nc'
V=squeeze(mean(double(ncread(pathdr,'V',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])),3));

Vdt=detrend(V',2)';  
Vp=V-Vdt;
stdVdt=std(Vdt,1,2);
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
eVme=stdVdt./sqrt(neff);





eV=squeeze(double(ncread(pathdr,'eV',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])));
%N=squeeze(double(ncread(pathdr,'eV',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])));


U=squeeze(mean(double(ncread(pathdr,'U',[stidxlat stidxlon 1],[lenidxlat lenidxlon 12])),3));
Udt=detrend(U',2)';  
Up=U-Udt;
stdUdt=std(Udt,1,2);
rr=[];
for i = 1:80
[r,l]=xcorr(Udt(i,:)','normalized');
rr(:,i)=r;
[m,lidx]=min(abs(r(284:304)));
lidx=lidx-1;
neffU(i)=284./lidx;
end
neffU=neffU';
neffU=smooth(neffU,5)';
neffU=smooth(neffU,5)';
eUme=stdUdt./sqrt(neffU);


lon=squeeze(double(ncread(pathdr,'Lon',[stidxlon],[lenidxlon])));
lat=squeeze(double(ncread(pathdr,'Lat',[stidxlat],[lenidxlat])));
Vme=mean(V,2,'omitnan');
Ume=mean(U,2,'omitnan');
for ij=1:80
    clear p
    p=polyfit(lon-mean(lon),U(ij,:),1);
    dUdxme(ij)=p(1)./110750;
end
dVdy=zeros(size(Vme));
dVdy(2:end-1)=(Vme(3:end)-Vme(1:end-2))./(0.5.*111555);


figure('position',[50 50 1200 400]);
subplot(1,2,1),...
%yyaxis left
plot(lat,Vme,'b-','linewidth',2);
hold on
ar=area(lat,[Vme-eVme 2.*eVme]);
ar(1).FaceAlpha=0;
ar(2).FaceAlpha=0.5;
ar(2).FaceColor='b'

grid on
xlim([-8 8])
ylim([-.1 .1])


zmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','depth'));






umit3d=squeeze((ncread('./data/mitgcm20yr_buoyavg.nc','u')));

dUdxmit=squeeze(mean(umit3d((end-20):end,:,:),1)-mean(umit3d(1:20,:,:),1))./(70.*110000);

umit=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','u'),1,'omitnan'));
umit2d=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','u'),1,'omitnan'));

wmit=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','w'),1,'omitnan'));
vmit=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','v'),1,'omitnan'));

vmit=mean(vmit(:,6:7),2,'omitnan');
umit=mean(umit(:,6:7),2,'omitnan');
wmit=mean(wmit(:,6:7),2,'omitnan');
latmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','latitude'));
lonmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','longitude'));



dVdymit=zeros(size(vmit));
dVdymit(2:end-1)=(vmit(3:end)-vmit(1:end-2))./(2.*0.0501251.*111555);
%yyaxis left
hold on;
plot(latmit,vmit,'r-','linewidth',2);
ylabel('m/s')
set(gca,'fontsize',13)
hold on;

set(gca,'fontsize',13)
set(gcf,'color','w')
title('Climatological meridional velocity at 15 m across the equator')


% Create textarrow
annotation(gcf,'textarrow',[0.362499999999999 0.426785714285714],...
    [0.380952380952382 0.32857142857143],...
    'Color','r',...
    'String',{'V at 15 m (MITgcm)'},'Fontsize',14);
% 
% Create textarrow
annotation(gcf,'textarrow',[0.582142857142857 0.510714285714286],...
    [0.299 0.302380952380952],'Color','b',...
    'String',{'V at 15 m (drifters)'},'Fontsize',14);


%print(gcf,'GDPclim_vs_MITgcmclim_V.png','-dpng','-r200')

%%
subplot(1,2,2),...
%yyaxis left
    plot(latmit(201:400),squeeze(dVdymit(201:400)).*86400,'b-','linewidth',3)
hold on;
    plot(-latmit(1:200),squeeze(dVdymit(1:200)).*86400,'r--','linewidth',3)
hold on;
    plot(latmit(201:400),squeeze(dVdymit(201:400)).*86400-...
        squeeze(dVdymit(200:-1:1)).*86400,'k:','linewidth',3)

    plot(lat(41:80),dVdy(41:80).*86400,'b-','linewidth',1)
hold on;
    plot(-lat(1:40),squeeze(dVdy(1:40)).*86400,'r--','linewidth',1)
hold on;
    plot(lat(41:80),squeeze(dVdy(41:80)).*86400-...
        squeeze(dVdy(40:-1:1)).*86400,'k:','linewidth',1)
xlabel('Latitude (degrees off equator)')
set(gcf,'color','w')
set(gca,'fontsize',15)

xlim([0.5 10])
ylim([-.02 .05])
grid on
ylabel('dV/dy (d^{-1})')
legend('North of the equator','South of the equator','Asymmetry at each latitude (North minus South)')
%yyaxis right
%ylim([-.02 .05].*25)
%ylabel('Implied vertical velocity at 50 m (m/d)')
title('Climatological divergence of meridional velocity dV/dy at 15 m')
%print(gcf,'GDPclim_vs_MITgcmclim_dVdy.png','-dpng','-r200')


%% maps of w50, dUdx15m+dVdy15m, dVdy15m
wmitmap=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','w')));
wmitmap=mean(wmitmap(:,:,20:21),3);
vmitmap=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','v')));
vmitmap=mean(vmitmap(:,:,6:7),3);
umitmap=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','u')));
umitmap=mean(umitmap(:,:,6:7),3);

load('./data/mitgcm_grid.mat','XC','YC');
load('./data/mitgcm_offline_KE_avg.mat');
V2=squeeze(mean(SFnow(:,:,6:7,2),3));
UV=squeeze(mean(SFnow(:,:,6:7,4),3));
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

dvdymap=zeros(size(vmitmap));
dvdymap(:,2:end-1)=(vmitmap(:,3:end)-vmitmap(:,1:end-2))./(2.*DYC(:,2:end-1));
dvdymap(:,1)=dvdymap(:,2);
dvdymap(:,end)=dvdymap(:,end-1);
dudxmap=zeros(size(vmitmap));
dudxmap(2:end-1,:)=(umitmap(3:end,:)-umitmap(1:end-2,:))./(2.*DXC(2:end-1,:));
dudxmap(1,:)=dudxmap(2,:);
dudxmap(end,:)=dudxmap(end-1,:);

dUVdy=zeros(size(UV));
dUVdy(:,2:end-1)=(UV(:,3:end)-UV(:,1:end-2))./(2.*DYC(:,2:end-1));
dUVdy(:,1)=dUVdy(:,2);
dUVdy(:,end)=dUVdy(:,end-1);








wmit=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','w'),1,'omitnan'));
Tmit=squeeze(mean(ncread('./data/mitgcm20yr_buoyavg.nc','theta'),1,'omitnan'));
zmit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','depth'));
ymit=squeeze(ncread('./data/mitgcm20yr_buoyavg.nc','latitude')).*111300;




wmit3d=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','w')));
wmit3dsm=zeros(size(wmit3d));
for j = 1:400
    j/400
    pause(0.01)
wmitmapsm(:,j)=smooth(wmitmap(:,j),141);
for k = 1:136
 wmit3dsm(:,j,k)=smooth(wmit3d(:,j,k),141);
end

end
%wmit=mean(wmit(:,20),2,'omitnan');
figure('position',[50 50 1200 1000]);
subplot(3,2,1),...
contourf(ymit./1000./111,zmit(2:end),86400.*wmit(:,2:end)',linspace(-2,2,21)); shading flat
hold on
caxis([-2 2])
caxis(caxis)
hold on
[c,h]=contour(ymit./1000./111,zmit,Tmit',-2:2:30,'linewidth',1,'color','k');
clabel(c,h);
ylabel('depth (m)','FontSize',14,'FontWeight','normal');
xlabel('latitude ({\circ}N)','FontSize',14,'FontWeight','normal');
cbh = colorbar();
set(gca,'FontSize',14,'FontWeight','normal');
title('(a) Zonal and time mean vertical velocity (m/d)','fontsize',14,'FontWeight','normal');
ylim([-300 0])
hold on;
plot(2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
plot(-2.*ones(4,1),linspace(-300,0,4),'k-','linewidth',2)
xlim([-6 6])
plot(linspace(-8,8,17),-50.*ones(17,1),'k-','linewidth',2)
grid on
set(gca,'fontsize',14)

colormap(cmocean('balance',20))
subplot(3,2,2),...
hold on;
plot(latmit(201:400),squeeze(mean(wmit(201:400,20:21),2)).*86400-...
    squeeze(mean(wmit(200:-1:1,20:21),2)).*86400,'k:','linewidth',2);
%plot(latmit(201:400),dMeydy(201:400).*86400-dMeydy(200:-1:1).*86400,'k-','linewidth',1)
hold on;
xlabel('latitude (degrees from equator)')
        plot(latmit(201:400),35.*(squeeze(dVdymit(201:400)).*86400-...
        squeeze(dVdymit(200:-1:1)).*86400),'-.','linewidth',2,'color','k')
    
plot(lat(41:80),35.*(squeeze(dVdy(41:80)).*86400-...
        squeeze(dVdy(40:-1:1)).*86400),'-.','linewidth',2,'Color',[.5 .5 .5])
plot(2.*ones(4,1),linspace(-0.05.*25,0.05.*25,4),'k-','linewidth',2)
plot(latmit(201:400),squeeze(mean(wmit(201:400,20:21),2)).*86400,'b-','linewidth',2)
hold on;
plot(-latmit(1:200),squeeze(mean(wmit(1:200,20:21),2)).*86400,'r-','linewidth',2)


xlim([0 6])
ylim([-.05 .05]*25)
%ylim([-.1 .1])
grid on
ylabel('(m/d)')
title('(b) Vertical velocity at 50 m by distance from the equator','fontweight','normal','fontsize',14)
%legend('North of the equator','South of the equator','Asymmetry by latitude (NH minus SH)','Asymmetry in Ekman suction')
legend('Asymmetry (North H. minus South H.)','Asymmetry H dV/dy','Observed asymmetry H dV/dy','location','southwest')

set(gcf,'color','w')
set(gca,'fontsize',14)

subplot(3,2,[3 4]),contourf(XC,YC,wmitmap.*86400,linspace(-4,4,41),'linestyle','none');
caxis([-2 2]);
colormap(cmocean('balance',20));
hold on
plot(XC(:,10),zeros(1420,1),'k-','linewidth',2)
plot(XC(:,10),2.*ones(1420,1),'k-','linewidth',2)
plot(XC(:,10),-2.*ones(1420,1),'k-','linewidth',2)
grid on;
ylim([-6 6])
cbh=colorbar;
ylabel(cbh,'m/d')
title('(c) Map of mean vertical velocity at 50 m','fontweight','normal','fontsize',14)
set(gca,'color','k','linewidth',1)
ylabel('latitude ({\circ}N)','FontSize',14,'FontWeight','normal');
xlabel('longitude ({\circ}E)','FontSize',14,'FontWeight','normal');
set(gca,'fontsize',14)
hold on
wmitmapsym(:,1:200)=wmitmap(:,1:200).*0.5+wmitmap(:,400:-1:201).*.5;
wmitmapsym(:,201:400)=wmitmap(:,200:-1:1).*0.5+wmitmap(:,201:400).*.5;
wmitmapasym=wmitmap-wmitmapsym;
wmitmapasym2(:,201:400)=wmitmapasym(:,201:400)-wmitmapasym(:,200:-1:1);
[w50latmax,idxmax]=max(wmitmap(:,101:300),[],2);
[w50asymlatmax,idxmax2]=max(wmitmapasym2(:,201:300),[],2);
lattemp=latmit(201:400);
plot(XC(1:5:end,200),latmit(100+idxmax(1:5:end)),'r.','markersize',5)
plot(XC(1:5:end,200),lattemp(idxmax2(1:5:end)),'b.','markersize',5)
plot(XC(1:5:end,200),-lattemp(idxmax2(1:5:end)),'b.','markersize',5)

set(gca,'XTick',[-160:10:-100])
set(gca,'XTickLabels',{'160','150','140','130','120','110','100'})
xlim([-168 -97])
xlabel('Longitude (^oW)')





[weqmax,idxeqmax]=max(mean(wmit3dsm(:,200:201,:),2),[],3);


wmitmapsym=zeros(size(wmitmap));
wmitmapsym(:,1:200)=wmitmapsm(:,1:200).*0.5+wmitmapsm(:,400:-1:201).*.5;
wmitmapsym(:,201:400)=wmitmapsm(:,200:-1:1).*0.5+wmitmapsm(:,201:400).*.5;
wmitmapasym=wmitmapsm-wmitmapsym;
wmitmapasym2(:,201:400)=wmitmapasym(:,201:400)-wmitmapasym(:,200:-1:1);
[w50latmax,idxmax]=max(wmitmapsm(:,101:300),[],2);
[w50asymlatmax,idxmax2]=max(wmitmapasym2(:,201:300),[],2);
lattemp=latmit(201:400);


%avg5Sto5N 50 m (lon)
subplot(3,2,[5 6],'Position',[0.13 0.11 0.743333333333333 0.215735294117647]),...
%plot(lonmit,86400.*smooth(sum(DXC(:,101:300).*wmitmap(:,101:300),2,'omitnan')./sum(DXC(:,101:300),2,'omitnan'),101)); hold on;
plot(lonmit,86400.*sum(DXC(:,161:240).*wmitmapsm(:,161:240),2,'omitnan')./sum(DXC(:,161:240),2,'omitnan'),'linewidth',1); hold on
ylabel('w (m/d)')
xlabel('Longitude')
%int2Sto2N 50 m (lon)

% w50max (lon)
plot(lonmit,86400.*w50latmax,'linewidth',1); 
plot(lonmit,86400.*w50asymlatmax,'linewidth',1); 
plot(lonmit,86400.*weqmax,'linewidth',1); 

grid on
title('(d) Zonal dependence of mean upwelling','fontweight','normal')
%legend('w_{50} averaged 2S-2N','max w_{50} at each lat','max asymmetry in w_{50} at each lat','max w on the equator at each lat','location','southwest')
% wmaxEq (lon)
set(gcf,'color','w')
set(gca,'fontsize',14)
ylim([0 2])
set(gca,'XTick',[-160:10:-100])
set(gca,'XTickLabels',{'160','150','140','130','120','110','100'})
xlim([-168 -97])
xlabel('Longitude (^oW)')

annotation(gcf,'textarrow',[0.484166666666667 0.496666666666667],...
    [0.446841303364005 0.532127017649719],'Color',[1 0 0],...
    'String',{'Latitude of Maximum Upwelling'},...
    'FontSize',18,...
    'FontName','Helvetica Neue');

% Create textarrow
annotation(gcf,'textarrow',[0.710000000000003 0.688333333333333],...
    [0.434613063810528 0.48363252375924],'Color',[0 0 1],...
    'String',{'Latitude of Maximum Asymmetry'},...
    'FontSize',18,...
    'FontName','Helvetica Neue');

% Create arrow
annotation(gcf,'arrow',[0.710833333333337 0.684166666666667],...
    [0.435663146779301 0.5480464625132],'Color',[0 0 1]);

% Create textarrow
annotation(gcf,'textarrow',[0.435833333333333 0.436666666666667],...
    [0.160506863780359 0.208025343189018],...
    'String',{'w_{50} averaged 2{\circ}S-2{\circ}N'},...
    'FontSize',14,...
    'FontName','Helvetica Neue');

% Create textarrow
annotation(gcf,'textarrow',[0.435 0.404166666666667],...
    [0.305174234424498 0.287222808870116],'String',{'max w_{50}'},'FontSize',14,...
    'FontName','Helvetica Neue');

% Create textarrow
annotation(gcf,'textarrow',[0.74952380952381 0.669166666666667],...
    [0.305117061396896 0.287222808870116],'String',{'max w on equator'},...
    'FontSize',14,...
    'FontName','Helvetica Neue');

% Create textarrow
annotation(gcf,'textarrow',[0.685833333333334 0.6875],...
    [0.142555438225977 0.178458289334741],'String',{'max asymmetry in w_{50}'},...
    'FontSize',14,...
    'FontName','Helvetica Neue');


print(gcf,'MITgcmclim_W_map_4panel.png','-dpng','-r200')








%%
Ujohnson=ncread('./data/johnson-eq-pac-adcp.nc','UM');

PDjohnson=ncread('./data/johnson-eq-pac-adcp.nc','SIGMAM');
XLONjohnson=ncread('./data/johnson-eq-pac-adcp.nc','XLON');
XLONjohnson=XLONjohnson(5:10);
YLATjohnson=ncread('./data/johnson-eq-pac-adcp.nc','YLAT11_101');
dUdxjohnson=zeros(91,50);
for i = 1:91
    for j=1:50
        clear p
        p=polyfit(squeeze(XLONjohnson-mean(XLONjohnson)),squeeze(Ujohnson(5:10,i,j)),1);
        dUdxjohnson(i,j)=p(1)./110000;
        clear p
        p=polyfit(XLONjohnson-mean(XLONjohnson),Ujohnson(5:10,i,j),3);
        Upjohnson(i,j)=mean(polyval(p,linspace(XLONjohnson(1)-mean(XLONjohnson),XLONjohnson(end)-mean(XLONjohnson),201)));
    end
end
ZDEPjohnson=ncread('./data/johnson-eq-pac-adcp.nc','ZDEP1_50');
Ujohnson=squeeze(mean(Ujohnson(5:10,:,:),1,'omitnan'));
PDjohnson=squeeze(mean(PDjohnson(5:10,:,:),1,'omitnan'));








%% T,S,rho structure
g = 9.81; %gravity
rhoref = 1035; % rho_0

LATargo=ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','LATITUDE',[331],[120]);
LONargo=ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','LONGITUDE',[1033],[426]);

Pargo=double(repmat(reshape(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','PRESSURE'),[1 1 58]),[426 120 1]));

Targo=double(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','ARGO_TEMPERATURE_MEAN',[1033 331 1],[426 120 58]));

Sargo=double(ncread('./data/RG_ArgoClim_33pfit_2019_mean.nc','ARGO_SALINITY_MEAN',[1033 331 1],[426 120 58]));

%
addpath ./seawater/
PDargo=sw_pden(Sargo,Targo,Pargo,zeros(size(Pargo)));
THargo=sw_ptmp(Sargo,Targo,Pargo,zeros(size(Pargo)));
Zargo=-sw_dpth(squeeze(Pargo(1,:,:)),repmat(LATargo,[1 58])); 
Svan=sw_svan(Targo, Sargo, Pargo);  
Dargo = calculate_dynamic_height_seawater(Targo, Sargo, squeeze(Pargo(1,1,:)));
Dargo500mref=-Dargo(:,:,1:34)+repmat(squeeze(Dargo(:,:,34)),[1 1 34]);
dargovelscale=9.81./(14.6e-5.*sind(-5))./(66.*110750);
figure;
subplot(2,1,1),...
contourf(repmat(LATargo,[1 34]),Zargo(:,1:34),dargovelscale.*squeeze(nanmean(Dargo500mref((end-29):end,:,:),1)-nanmean(Dargo500mref(1:30,:,:),1)),-2:.5:10); 
xlim([-10 10])
caxis([-1 6])
ylim([-300 0])
colormap(jet)
cbh=colorbar;
ylabel(cbh,'cm/s')
xlabel('latitude')
ylabel('depth (m)')
title('(a) Zonal dynamic height difference scaled by g/fL at 5^{\circ}N (Argo)')
Tmit3d=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','theta')));
Smit3d=squeeze(double(ncread('./data/mitgcm20yr_buoyavg.nc','salt')));
PDmit=sw_dens(Smit3d,Tmit3d,zeros(size(Smit3d)));




PDmit2d=squeeze(mean(PDmit,1,'omitnan'));
Smit2d=squeeze(mean(Smit3d,1,'omitnan'));
Tmit2d=squeeze(mean(Tmit3d,1,'omitnan'));

THargo2d=squeeze(mean(THargo,1,'omitnan'));
Sargo2d=squeeze(mean(Sargo,1,'omitnan'));
PDargo2d=squeeze(mean(PDargo,1,'omitnan'));

dmitvelscale=9.81./(14.6e-5.*sind(-5))./(66.*110750);

Pmit=sw_pres(abs(zmit),0);
ITmit3d = sw_temp(Smit3d,Tmit3d,repmat(reshape(Pmit,[1 1 136]),[1420 400 1]),zeros(size(Tmit3d)));
Dmit3d = calculate_dynamic_height_seawater(ITmit3d, Smit3d, squeeze(Pmit));
clear ITmit3d
Dmit500mref=-Dmit3d(:,:,1:125)+repmat(squeeze(mean(Dmit3d(:,:,125:126),3)),[1 1 125]);

 

subplot(2,1,2),...
contourf(latmit,zmit(1:125),dmitvelscale.*squeeze(nanmean(Dmit500mref((end-100):end,:,:),1)-nanmean(Dmit500mref(1:100,:,:),1))',-2:.5:10); 
xlim([-10 10])
caxis([-1 6])
colormap(jet)
cbh=colorbar;
ylim([-300 0])
title('(b) MITgcm')
ylabel(cbh,'cm/s')
xlabel('latitude')
ylabel('depth (m)')
set(gcf,'color','w')
%%
figure('position',[50 50 1300 800]);
subplot(2,3,1),...
contourf(lonmit,latmit,mean(Tmit3d(:,:,1:2),3)',20:.5:32);  cbh=colorbar; ylabel(cbh,'deg C');
caxis([24 30])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean T at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(12))

subplot(2,3,4),...
contourf(LONargo-360,LATargo,mean(THargo(:,:,1),3)',20:.5:32);  cbh=colorbar; ylabel(cbh,'deg C');
caxis([24 30])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean S at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(12))

subplot(2,3,2),...
contourf(lonmit,latmit,mean(Smit3d(:,:,1:2),3)',33:.1:36);  cbh=colorbar; ylabel(cbh,'psu');
caxis([33 36])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean S at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(30))

subplot(2,3,5),...
contourf(LONargo-360,LATargo,mean(Sargo(:,:,1),3)',33:.1:36);  cbh=colorbar; ylabel(cbh,'psu');
caxis([33 36])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean S at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(30))

subplot(2,3,3),...
contourf(lonmit,latmit,mean(PDmit(:,:,1:2)-1000,3)',20:.25:30);  cbh=colorbar; ylabel(cbh,'kg/m3');
caxis([21 25])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean PD at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(16))

subplot(2,3,6),...
contourf(LONargo-360,LATargo,mean(PDargo(:,:,1)-1000,3)',20:.25:30);  cbh=colorbar; ylabel(cbh,'kg/m3');
caxis([21 25])
ylim([-10 10])
xlim([-167.5 -97.5])
title('MITgcm mean PD at 2.5 m')
xlabel('longitude')
ylabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(16))
set(gcf,'color','w')

print(gcf,'SSTSSSSPD_map_mitgcmvsArgo.png','-dpng','-r200')


figure('position',[50 50 1300 800]);
subplot(2,3,1),...
contourf(latmit,zmit(1:112),Tmit2d(:,1:112)',10:1:32);  cbh=colorbar; ylabel(cbh,'deg C');
caxis([10 30])
xlim([-10 10])
ylim([-300 0])
title('(a) MITgcm zonal mean T')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(2,3,4),...
contourf(LATargo,mean(Zargo(:,1:26),1),THargo2d(:,1:26)',10:1:32);  cbh=colorbar; ylabel(cbh,'deg C');
caxis([10 30])
xlim([-10 10])
ylim([-300 0])
title('(b) Argo zonal mean T')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(2,3,2),...
contourf(latmit,zmit(1:112),Smit2d(:,1:112)',32:.1:36);  cbh=colorbar; ylabel(cbh,'psu');
caxis([34 36])
xlim([-10 10])
ylim([-300 0])
title('(c) MITgcm zonal mean S')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(2,3,5),...
contourf(LATargo,mean(Zargo(:,1:26),1),Sargo2d(:,1:26)',32:.1:36);  cbh=colorbar; ylabel(cbh,'psu');
caxis([34 36])
xlim([-10 10])
ylim([-300 0])
title('(d) Argo zonal mean S')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(2,3,3),...
contourf(latmit,zmit(1:112),-1000+PDmit2d(:,1:112)',20:.25:30);  cbh=colorbar; ylabel(cbh,'kg/m3');
caxis([21 27])
xlim([-10 10])
ylim([-300 0])
title('(e) MITgcm zonal mean PD')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(24))

subplot(2,3,6),...
contourf(LATargo,mean(Zargo(:,1:26),1),-1000+PDargo2d(:,1:26)',20:.25:30);  cbh=colorbar; ylabel(cbh,'kg/m3');
caxis([21 27])
xlim([-10 10])
ylim([-300 0])
title('(f) Argo zonal mean PD')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(24))

set(gcf,'color','w')

print(gcf,'THSPD_section_mitgcmvsArgo.png','-dpng','-r200')








for i = 1:120
    for j = 1:34
        clear p
        p=polyfit(squeeze(LONargo-mean(LONargo)),squeeze(Dargo500mref(:,i,j)),3);
        Dargo500mrefp(:,i,j)=polyval(p,LONargo-mean(LONargo));        
    end
end

[latga,longa]=meshgrid(LATargo,LONargo);
fga=7.29e-5.*2.*sind(latga);

DYCa=zeros(size(latga));
DYCa(:,2:end-1)=111300.*(latga(:,3:end)-latga(:,1:end-2))./2;
DYCa(:,1)=DYCa(:,2);
DYCa(:,end)=DYCa(:,end-1);

Ugargo = zeros(size(Dargo500mref));
Ugargo(:,2:end-1,:)=-g.*(Dargo500mrefp(:,3:end,:)-Dargo500mrefp(:,1:end-2,:))./(2.*DYCa(:,2:end-1))./fga(:,2:end-1)./1e2;
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

for j = 1:34
    clear p
    if ~isnan(utao140argo(j))
    p=polyfit([190 220 250]-mean(LONargo),[utao170argo(j) utao140argo(j) utao110argo(j)],2);
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


figure('position',[10 10 1100 600]);
subplot(3,2,1),...
contourf(latmit,zmit(1:112),umit2d(:,1:112)',-1:.1:1);  cbh=colorbar; ylabel(cbh,'ms^{-1}');
caxis([-1 1])
xlim([-8 8])
ylim([-300 0])
title('(a) MITgcm zonal mean U')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(3,2,2),...
contourf(latmit,zmit(1:112),1e8.*dUdxmit(:,1:112)',-20:1:20);  cbh=colorbar; ylabel(cbh,'10^{-8} s^{-1}');
caxis([-10 10])
xlim([-8 8])
ylim([-300 0])
title('(b) MITgcm zonal mean dU/dx')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))



subplot(3,2,3),...
contourf(YLATjohnson,-ZDEPjohnson(1:30),Upjohnson(:,1:30)',-1:.1:1);  cbh=colorbar; ylabel(cbh,'ms^{-1}');
caxis([-1 1])
xlim([-8 8])
ylim([-300 0])
title('(c) Johnson zonal mean U')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))

subplot(3,2,4),...
contourf(YLATjohnson,-ZDEPjohnson(1:30),1e8.*dUdxjohnson(:,1:30)',-20:1:20);  cbh=colorbar; ylabel(cbh,'10^{-8} s^{-1}');
caxis([-10 10])
xlim([-8 8])
ylim([-300 0])
title('(d) Johnson zonal mean dU/dx')
ylabel('depth (m)')
xlabel('latitude')
set(gca,'fontsize',13)
colormap(gca,jet(20))


subplot(3,2,5),...
contourf(LATargo',Zargo(10,1:34)',squeeze(mean(Ugargo,1))',linspace(-2,2,41));
colormap(jet(20));
caxis([-1 1]);
ylim([-300 0])
xlim([-8 8])
cbh=colorbar;
xlabel('latitude')
ylabel('depth (m)')
set(gcf,'color','w')
title('(e) Argo zonal mean U')
ylabel(cbh,'ms^{-1} ')
set(gca,'fontsize',14)

subplot(3,2,6),...
contourf(LATargo',Zargo(10,1:34)',(dUdxargo.*1e8)',linspace(-20,20,41));
colormap(jet(20));
caxis([-10 10]);
ylim([-300 0])
xlim([-8 8])
cbh=colorbar;
xlabel('latitude')
ylabel('depth (m)')
set(gcf,'color','w')
title('(f) Argo zonal mean dU/dx')
ylabel(cbh,'10^{-8} s^{-1} ')
set(gca,'fontsize',14)

set(gcf,'color','w')

print(gcf,'U_sections_mitgcmvsjohnsonvsargo.png','-dpng','-r200')


figure('position',[10 10 500 700]);
subplot(2,1,1),...
plot(lat,Ume,'k-','linewidth',1);
hold on
ar2=area(lat,[Ume-eUme 2.*eUme]);
ar2(1).FaceAlpha=0;
ar2(2).FaceAlpha=0.5;
ar2(2).FaceColor='k';
hold on;
plot(latmit,umit,'b--','linewidth',1);
hold on;
plot(YLATjohnson,Ujohnson(:,2),'r:','linewidth',1);
plot(LATargo,squeeze(mean(Ugargo(:,:,1),1)),'magenta-.','linewidth',1);
ylabel('m/s')
set(gca,'fontsize',13)
set(gcf,'color','w')
grid on
xlim([-8 8])
ylim([-.6 .6])
hold on
title('(a) Climatological zonal mean zonal velocity U at 15 m depth')
xlabel('latitude')

subplot(2,1,2),...
plot(lat,1e8.*dUdxme,'k-','linewidth',1);
hold on
plot(latmit,1e8.*squeeze(mean(dUdxmit(:,6:7),2)),'b--','linewidth',1);
hold on;
plot(YLATjohnson,1e8.*dUdxjohnson(:,2),'r:','linewidth',1);
plot(LATargo,1e8.*squeeze(mean(dUdxargo(:,2:3),2)),'magenta-.','linewidth',1);
ylabel('10^{-8} s^{-1}')
set(gca,'fontsize',13)
set(gcf,'color','w')
grid on
xlim([-8 8])
ylim([-8 8])
hold on
legend('Laurindo Gridded Drifters','MITgcm simulation','Johnson SADCP','Argo+MADCP')
xlabel('latitude')

title('(b) Climatological zonal mean dU/dx at 15 m depth')

print(gcf,'GDPclim_vs_MITgcmclim_UdUdx.png','-dpng','-r200')





