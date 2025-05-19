
%% Dan Whitt (daniel.b.whitt@nasa.gov) 
%% Copyright Dan Whitt
%% Written With Matlab v2023b 
%% Github upload May 19 2025
% dependencies
%~/'OneDrive - NASA'/documents/mitgcm_offline_KE_avg.mat
%~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc
%~/Downloads/adcp0n140w_dy.nc
clear all;
close all;
load ~/'OneDrive - NASA'/documents/mitgcm_offline_KE_avg.mat
zmit=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','depth'));
latmit=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','latitude'));
lonmit=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','longitude'));
umit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','u'));
vmit3d=squeeze(ncread('~/OneDrive - NASA/documents/mitgcm20yr_buoyavg.nc','v'));
umit0n140w=squeeze(umit3d(561,201,:));
vmit0n140w=squeeze(vmit3d(561,201,:));
clear umit3d vmit3d

lonmit(561)
latmit(201)
figure;
subplot(1,3,1),...
plot(squeeze(SFnow(561,201,:,1))-umit0n140w.^2,zmit);
ylim([-300 0])

subplot(1,3,2),...
plot(squeeze(SFnow(561,201,:,2))-vmit0n140w.^2,zmit);
ylim([-300 0])

subplot(1,3,3),...
plot(squeeze(SFnow(561,201,:,4))-umit0n140w.*vmit0n140w,zmit);
ylim([-300 0])

uobs=squeeze(double(ncread('~/Downloads/adcp0n140w_dy.nc','u_1205')));
vobs=squeeze(double(ncread('~/Downloads/adcp0n140w_dy.nc','v_1206')));
depth=double(ncread('~/Downloads/adcp0n140w_dy.nc','depth'));

U2obs=nanmean((uobs./100).^2,2);
Uobsm=nanmean(uobs./100,2);
V2obs=nanmean((vobs./100).^2,2);
Vobsm=nanmean(vobs./100,2);

UVobs=nanmean(vobs.*uobs./1e4,2);

subplot(1,3,1),...
    hold on;
plot(U2obs-Uobsm.^2,-depth);
ylim([-200 -30])
title('(a) <u^2>-<u>^2 m^2/s^2')
grid on
ylabel('Depth m')
set(gca,'fontsize',13)
xlim([0 1.5])

subplot(1,3,2),...
        hold on;

plot(V2obs-Vobsm.^2,-depth);
ylim([-200 -30])
title('(b) <v^2>-<v>^2 m^2/s^2')
grid on
set(gca,'fontsize',13)
xlim([0 .1])


subplot(1,3,3),...
        hold on;

plot(UVobs-Uobsm.*Vobsm,-depth);
ylim([-200 -30])
title('(c) <uv>-<u><v> m^2/s^2')
set(gca,'fontsize',13)
grid on
set(gcf,'color','w')
xlim([-.03 .03])

