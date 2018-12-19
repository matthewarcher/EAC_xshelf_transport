% Create Hovmoller of AVISO 2D Gridded SLA and U/V

% STEPS
% (1) Rotate u/v to the general coastline orientation
% (2) Isolate an angled line to pass a chosen longitude (custom choice)
% (3) Plot Hovmoller, plot time mean

% Written by: Matt Archer, early 2018

% LOAD DATA & CONCATENATE
load .\Post_Cruise\AVISO_1993_2005.mat; 
sla1 = sla;ug1=ug;vg1=vg;t1=t;clear ug vg t sla
load .\Post_Cruise\AVISO_2006_2018.mat;
sla2 = sla;ug2=ug;vg2=vg;t2=t;clear ug vg t sla
sla=cat(3,sla1,sla2);t =cat(1,t1,t2);ug=cat(3,ug1,ug2);vg=cat(3,vg1,vg2);
clear sla1 sla2 vg1 vg2 ug1 ug2 t1 t2 sl time_g
% Create 2D grid of longitude/latitude
LON = repmat(lon_sel,1,length(lat_sel));
LAT = repmat(lat_sel,1,length(lon_sel))';

% % Get Angle off North (for general coastline orientation)
% [east,north]=lonlat2km(LON(io2,ia2),LAT(io2,ia2),LON(io1,ia1),LAT(io1,ia1))
% tilt = rad2deg(atan(east/north));
tilt=27;

%% Rotate U and V to this angle
rot_deg_angle = -tilt; % for right oriention angle
uc = cosd(rot_deg_angle).*ug + sind(rot_deg_angle).*vg;  %  across-shelf
vc = -sind(rot_deg_angle).*ug + cosd(rot_deg_angle).*vg; %  along-shelf

%% Plot Figure

% Mean Sea Level Anomaly (sla)
SAM = nanmean(sla,3); % just to provide a visual

figure
% Set Axes etc.
axis square;set(gca,'fontsize',20)
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
box on
% Title and Labels
xlabel('Longitude','Fontsize',14), ylabel('Latitude','Fontsize',14), set(gca,'Fontsize',12)
set(gca,'fontsize',20)
ylim([-38 -30])
xlim([148 156])
hold on
pcolor(LON,LAT,SAM);shading flat
plot(LON,LAT,'k.')
hold on 

% Guiding Dots
[~,io1] = min(abs((lon_sel - 154.4)))%3)))
[~,io2] = min(abs((lon_sel - 150.6)))
[~,ia1] = min(abs((lat_sel - -30)))%2.57)))
[~,ia2] = min(abs((lat_sel - -36.4)))
hold on
plot(LON(io2,ia2),LAT(io2,ia2),'or','markerfacecolor','r');hold on
plot(LON(io1,ia1),LAT(io1,ia1),'or','markerfacecolor','r')

% Custom Choose the Line for Hovmoller (or click red dots)
[XP,YP] = ginput(2);
[~,iXPn] = min(abs(lon_sel-max(XP)))
[~,iYPn] = min(abs(lat_sel-max(YP)))
[~,iXPs] = min(abs(lon_sel-min(XP)))
[~,iYPs] = min(abs(lat_sel-min(YP)))

%% Get Data Along Chosen Transect

lati = find(lat_sel > lat_sel(iYPs) & lat_sel < lat_sel(iYPn))

lon_all = interp1([lat_sel(iYPs) lat_sel(iYPn)],...
    [lon_sel(iXPs) lon_sel(iXPn)],lat_sel(lati))

hold on
plot(lon_all,lat_sel(lati))

% Get longitude index
for i = 1:length(lati)
    [~,loni(i)]=min(abs(lon_sel-lon_all(i))); % I.e. get closest points
end

hold on
plot(lon_sel(loni),lat_sel(lati),'o','markersize',10,'markerfacecolor','w')

%% Produce the Hovmoller!

% Generate matrix from vector
tmat = repmat(t,1,length(lati))';
lmat = repmat(lat_sel(lati),1,length(t));

% Grab data
for k = 1:length(loni)  
    ucmat(k,1:length(t)) = squeeze(uc(loni(k),lati(k),:));   
    vcmat(k,1:length(t)) = squeeze(vc(loni(k),lati(k),:));
    slamat(k,1:length(t)) = squeeze(sla(loni(k),lati(k),:));
end

% Plot Sea Level Anomaly
figure;
subplot(2,1,1)
pcolor(tmat,lmat,slamat);shading flat
caxis([-.5 .5])
datetick
axis tight
set(gca,'fontsize',20)
ylabel('SLA (m)'), colorbar

% Plot Cross-Shelf Velocity
subplot(2,1,2)
pcolor(tmat,lmat,ucmat);shading flat
% hold on
% contour(tmat,lmat,vcmat,[0 .5 1],'k','linewidth',1);shading flat
datetick
caxis([-.5 .5])
rwbcmap = redblue(64);
colormap(rwbcmap)
axis tight
set(gca,'fontsize',20)
ylabel('Cross-Shelf Velocity (m s^-^1)'), colorbar

%% Take a Look at a Chosen Time

% Choose a time period example
ind = find(t == datenum(2017,09,12,0,0,0))
k=ind;

figure
pcolor(LON,LAT,sla(:,:,k));shading flat
hold on
quiver(LON,LAT,ug(:,:,k),vg(:,:,k),'k','linewidth',2)
ylim([-38 -30])
xlim([148 156])
title(datestr(t(k)))
set(gca,'fontsize',20)

hold on
hold on
plot(lon_sel(loni),lat_sel(lati),'o','markersize',10,'markerfacecolor','w')


%% Look at dipole effect on the mean fields in time

figure;plot(nanmean(ucmat,2)...
    -nanmean(nanmean(ucmat,2),1),lmat(:,1),'linewidth',3);set(gca,'fontsize',20)

figure;plot(nanmean(slamat,2),lmat(:,1),'linewidth',3);set(gca,'fontsize',20)

